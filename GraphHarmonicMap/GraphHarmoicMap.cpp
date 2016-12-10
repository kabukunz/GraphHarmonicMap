#include <cstdlib>
#include "GraphHarmoicMap.h"
#include <lemon/dijkstra.h>
#include <omp.h>
#include <queue>
#include <set>
#include "Parser/parser.h"

CGraphHarmoicMap::CGraphHarmoicMap()
{
    mesh = new CMesh();
    graph = new CGraph();
}


CGraphHarmoicMap::~CGraphHarmoicMap()
{
}

int CGraphHarmoicMap::setMesh(string filename)
{
    mesh->read_m(filename.c_str());
    return 0;
}

/*
* read a target graph, place source vertex randomly on the graph: arc is chosen randomly, always place target at the middle of each arc
*/
int CGraphHarmoicMap::setGraph(string filename)
{
    graph->read(filename);
    if (!mesh)
    {
        cerr << "read mesh first" << endl;
		exit(-1);
        return -1;
    }

	vector<SmartGraph::Edge> edges;
	for (SmartGraph::EdgeIt e(graph->g); e != INVALID; ++e)
	{
		edges.push_back(e);
	}

	for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
	{
		CVertex * v = *vit;
		string s = v->string();
		CParser parser(s);
		list<CToken*> & tokens = parser.tokens();
		
		for (auto tit = tokens.begin(); tit != tokens.end(); ++tit)
		{
			CToken * pt = *tit;
			if (pt->m_key == "target")
			{
				string str = strutil::trim(pt->m_value, "()");
				istringstream iss(str);
				int edgeId = 0;
				int nodeId = 0;
				double length = 0.0;
				iss >> edgeId >> nodeId >> length;				
				CTarget * t = new CTarget();
				t->edge = graph->g.edgeFromId(edgeId);				
				t->node = graph->g.u(t->edge);
				t->length = length;
				v->prop("target") = t;
			}
		}
	}
    
	return 0;

    int ne = edges.size();
    for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
    {
        CVertex * v = *vit;
        CTarget * t = new CTarget();
		double z = v->point()[2];
		if (z > 0.4)
		{
			t->edge = edges[2];
		}
		else if (z < -0.4)
		{
			t->edge = edges[1];
		} 
		else
		{
			t->edge = edges[0];
		}
		auto u = graph->g.u(t->edge);
		t->node = u;
		t->length = graph->edgeLength[t->edge] / 2;
        /*int i = rand() % ne;
        t->edge = edges[i];
		auto u = graph->g.u(t->edge);
        t->node = u;
		double r = rand() / (double)RAND_MAX;
        t->length = graph->edgeLength[t->edge] * r;*/
        v->prop("target") = t;
    }
    return 0;
}

int CGraphHarmoicMap::calculateEdgeLength()
{
    for (MeshEdgeIterator eit(mesh); !eit.end(); ++eit)
    {
        CEdge * e = *eit;
        CVertex * v1 = e->halfedge(0)->source();
        CVertex * v2 = e->halfedge(0)->target();
        e->length() = (v1->point() - v2->point()).norm();
    }
    return 0;
}

int CGraphHarmoicMap::calculateEdgeWeight()
{
    auto inverse_cosine_law = [](double li, double lj, double lk) { return acos((lj * lj + lk * lk - li * li) / (2 * lj * lk)); };
    for (MeshEdgeIterator eit(mesh); !eit.end(); ++eit)
    {
        CEdge * e = *eit;
        e->prop("weight") = double(0.0);
        e->halfedge(0)->touched() = false;
        e->halfedge(1)->touched() = false;
    }
    for (MeshHalfEdgeIterator heit(mesh); !heit.end(); ++heit)
    {
        CHalfEdge * he = *heit;
        if (he->touched()) continue;
        CHalfEdge * he_next = he->he_next();
        CHalfEdge * he_prev = he->he_prev();
        double theta1 = inverse_cosine_law(he->edge()->length(), he_next->edge()->length(), he_prev->edge()->length());
        he->edge()->prop("weight") = 0.5 / tan(theta1);
        he->touched() = true;

        CHalfEdge * he_dual = he->he_sym();
        if (!he_dual || he_dual->touched()) continue;
        CHalfEdge * he_dual_next = he_dual->he_next();
        CHalfEdge * he_dual_prev = he_dual->he_prev();
        double theta2 = inverse_cosine_law(he_dual->edge()->length(), he_dual_next->edge()->length(), he_dual_prev->edge()->length());
        he->edge()->prop("weight") = 0.5 / tan(theta1) + 0.5 / tan(theta2);
        he_dual->touched() = true;
    }
    return 0;
}

/*
 * find the minimum point of a quadratic function in interval [x0, x1]
 */
inline double quadraticMininum(double a, double b, double c, double x0, double x1, double &x)
{
    auto func = [=](double x) { return a * x * x + b * x + c; };
	/*if (fabs(x0 - x1) < EPS)
	{
		x = (x0 + x1) / 2;
		return func(x);
	}*/
    double f0 = func(x0);
    double f1 = func(x1);
    x = x0;
    double fm = f0;
    if(f1 < f0)
    {
        x = x1;
        fm = f1;
    }
    double x2 = -b/a/2.0;
    double f2 = func(x2);
    if(x0 < x2 && x2 < x1)
    {
        if(f2 < fm)
        {
            x = x2;
            fm = f2;
        }
    }
    return fm;
}

/*
 * distance between two points x and y on the graph 
 */
double CGraphHarmoicMap::distance(CTarget * x, CTarget * y)
{
    // first compute the shorest distance between the node of two edges
    auto & ex = x->edge;
    auto & ey = y->edge;
	double elx = graph->edgeLength[ex];
	double ely = graph->edgeLength[ey];
	auto & ux = graph->g.u(ex);
	auto & vx = graph->g.v(ex);
	auto & uy = graph->g.u(ey);
	auto & vy = graph->g.v(ey);
	SmartGraph::Node nx, ny;
	double dey = distance(x, ey, nx, ny);
	if (uy == vy) // y is a loop
	{
		if (ex == ey)
		{
			dey = fabs(x->length - y->length);
			if (dey > ely / 2.0) dey = ely - dey;
		} 
		else
		{
			if (y->length > ely / 2.0) dey += ely - y->length;
			else dey += y->length;
		}
	}
	else
	{
		if (ex == ey) // x, y on same edge
		{
			dey = fabs(x->length - y->length);
		} 
		else
		{
			// ny is the node nearer to x, if it is the starting node of y
			if (ny == uy) dey = fabs(dey - y->length);
			else dey += ely - y->length;
		}
		
	}
	return dey;
}

double CGraphHarmoicMap::distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne)
{
	auto u = graph->g.u(e);
	auto v = graph->g.v(e);

	SmartGraph::Node ux, vx;
	if (u == v) // e is a loop
	{
		double du = distance(x, u, ux);
		nx = ux;
		ne = u;
		return du;
	}
	
	double du = distance(x, u, ux);
	double dv = distance(x, v, vx);
	if (du < dv)
	{
		nx = ux;
		ne = u;
		return du;
	}
	else
	{
		nx = vx;
		ne = v;
		return dv;
	}
}

double CGraphHarmoicMap::distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx)
{
	bool isLoop = false;
	auto e = x->edge;
	double el = graph->edgeLength[e];	

	auto u = graph->g.u(e);
	auto v = graph->g.v(e);
		
	SmartGraph::NodeMap<double> dist(graph->g);
	auto dij = dijkstra(graph->g, graph->edgeLength).distMap(dist);
	
	if (u == v)
	{
		double dx = x->length;
		if (dx > el / 2.0) dx = el - dx;
		double d = 0.0;
		dij.dist(d).run(n, u);
		nx = u;
		return dx + d;
	}
	else 
	{
		double d1 = 0.0, d2 = 0.0;
		dij.dist(d1).run(n, u);
		dij.dist(d2).run(n, v);
		if (d1 < d2)
		{
			double dx = x->length;
		 	nx = u;
			return dx + d1;
		}
		else
		{
			double dx = el - x->length;
			nx = v;
			return dx + d2;
		}
	}	
}

double CGraphHarmoicMap::distance(CTarget * x, SmartGraph::Node n)
{
	SmartGraph::Node nx;
	return distance(x, n, nx);
}


double CGraphHarmoicMap::calculateBarycenter(CVertex * v, vector<CVertex*> & nei)
{
    //vector<CVertex*> nei;	
	vector<CTarget*> neit;
    vector<double> ew;
	for (auto * vj : nei)
	{
		void * t = vj->prop("target");
		CTarget * vt = (CTarget*)t;
		neit.push_back(vt);
		double w = mesh->vertexEdge(v,vj)->prop("weight");
		ew.push_back(w);
	}
	/*for (VertexOutHalfedgeIterator heit(mesh, v); !heit.end(); ++heit)
    {
        CHalfEdge * he = *heit;
		CVertex * vj = he->target();
		void * t = vj->prop("target");
		CTarget * vt = (CTarget*)t;
        nei.push_back(vj);
		neit.push_back(vt);
		double w = he->edge()->prop("weight");
        ew.push_back(w);
    }*/
    double a = 0.0;
    for (double w : ew) a += w;
    
    vector<double> fm;
    vector<CTarget*> vx;
    for (SmartGraph::EdgeIt e(graph->g); e != INVALID; ++e)
    {		
		auto u = graph->g.u(e);
		auto v = graph->g.v(e);
		double el = graph->edgeLength[e];
		vector<bool> be;
		vector<double> bp = { 0, el / 2.0, el };
		for (auto it = neit.begin(); it != neit.end(); ++it)
		{
			CTarget * vt = *it;
			bool b = vt->edge == e;
			if (b)
			{
				if (vt->length < el / 2.0) bp.push_back(vt->length + el / 2.0);
				else bp.push_back(vt->length - el / 2.0);
			}
			be.push_back(b);
		}
        if (u == v) // if e is a loop
        {
			sort(bp.begin(), bp.end());
			vector<double> bx0;
			for (int i = 0; i < neit.size(); ++i)
			{
				if (be[i])
				{
					double x = -neit[i]->length;
					bx0.push_back(x);
				}
				else
				{
					double dy = distance(neit[i], u);
					bx0.push_back(dy);
				}
			}
			// find minimum point in this loop
			CTarget * t = new CTarget();
			t->edge = e;
			t->node = u;
			t->length = 0.0;
			double ei = 1.0 / EPS;
			for (int i = 0; i < bp.size() - 1; ++i)
			{
				vector<double> bx(bx0);
				double x0 = bp[i];
				double x1 = bp[i + 1];
				double xm = (x0 + x1) / 2.0;
				
				double b = 0.0, c = 0.0;
				for (int j = 0; j < bx.size(); ++j)
				{
					if (be[j])
					{
						if (xm + bx[j] > el / 2.0) bx[j] -= el;
						else if (xm + bx[j] < -el / 2.0) bx[j] += el;
					}
					else
					{
						if (xm > el / 2.0) bx[j] = -bx[j] - el;
					}
					b += 2 * ew[j] * bx[j];
					c += ew[j] * bx[j] * bx[j];
				}
				double x = 0;
				double mi = quadraticMininum(a, b, c, x0, x1, x);
				if (mi < ei)
				{
					ei = mi;
					t->length = x;
				}
			}			
			fm.push_back(ei);
			vx.push_back(t);            
        }
        
        else // if e is not loop
        {
			vector<double> bx;
			for (int i = 0; i < neit.size(); ++i)
			{
				if (be[i])
				{
					double x = -neit[i]->length;
					bx.push_back(x);
				}
				else
				{
					double du = distance(neit[i], u);
					double dv = distance(neit[i], v);
					if (du <= dv) bx.push_back(du);
					else bx.push_back(-dv - el);
				}
			}
			double b = 0.0, c = 0.0;
			for (int j = 0; j < bx.size(); ++j)
			{
				b += 2 * ew[j] * bx[j];
				c += ew[j] * bx[j] * bx[j];
			}
			CTarget * t = new CTarget();
			t->edge = e;
			t->node = u;
			t->length = 0.0;
			double x = 0;
			double mi = quadraticMininum(a, b, c, 0, el, x);
			t->length = x;
			fm.push_back(mi);
			vx.push_back(t);
        }
        
    }

	double fmm = fm[0];
	CTarget * vm = vx[0];
	for (int i = 1; i < fm.size(); ++i)
	{
		if (fm[i] < fmm)
		{
			fmm = fm[i];
			vm = vx[i];
		}
	}
	void * t = v->prop("target");
	CTarget * vt = (CTarget*)t;
	double dv = distance(vm, vt);
	vt->edge = vm->edge;
	vt->node = vm->node;
	vt->length = vm->length;
	
	for (int i = 0; i < vx.size(); ++i)
	{
		delete vx[i];
	}
    return dv;
}

/*
 * for every vertex, move it to its neighbors' weighted center
*/
int CGraphHarmoicMap::harmonicMap()
{
    calculateEdgeLength();
    calculateEdgeWeight();
	vector<CVertex*> vv;
	vector<vector<CVertex*>> neis;
	
	for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
	{
		CVertex * v = *vit;
		vv.push_back(v);
		vector<CVertex*> nei;
		for (VertexOutHalfedgeIterator heit(mesh, v); !heit.end(); ++heit)
		{
			CHalfEdge * he = *heit;
			CVertex * vj = he->target();						
			nei.push_back(vj);
		}
		neis.push_back(nei);
	}

	int k = 0;
    while (k++ < 5000)
    {
		double err = 0;
		//random_shuffle(vv.begin(), vv.end());
		//#pragma omp parallel for
		for(int i = 0; i < vv.size(); ++i)
        {		
            double d = calculateBarycenter(vv[i], neis[i]);
			if (d > err) err = d;
        }
		if (k % 100 == 0) cout << "#" << k << ": " << err << endl;
		if (k % 500 == 0) writeMap("harmonic." + to_string(int(k / 500)));
		if (err < EPS) break;
    }

    return 0;
}

void CGraphHarmoicMap::test()
{
	auto v0 = mesh->idVertex(8);
	auto v1 = mesh->idVertex(100);
	void * x = v0->prop("target");	
	void * y = v1->prop("target");
	CTarget * tx = (CTarget*)x;
	CTarget * ty = (CTarget*)y;
	double d = distance(tx, ty);

	calculateEdgeLength();
	calculateEdgeWeight();
	//double db = calculateBarycenter(v1);
	harmonicMap();
}

int CGraphHarmoicMap::writeMap(string filename)
{
	ofstream map(filename);
	for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
	{
		CVertex * v = *vit;
		void * t = v->prop("target");
		CTarget * vt = (CTarget*)t;
		CPoint & p = v->point();
		map << graph->g.id(vt->edge) << " " << graph->g.id(vt->node) << " " << vt->length << endl;
	}

    return 0;
}
