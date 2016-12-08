#include "GraphHarmoicMap.h"



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

    int ne = edges.size();
    for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
    {
        CVertex * v = *vit;
        CTarget * t = new CTarget();
        int i = rand() % ne;
        t->edge = edges[i];
        t->direction = rand() % 2;
        t->length = graph->edgeLength[t->edge] / 3.0;
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
    auto cosine_law = [](double li, double lj, double lk) { return acos((lj * lj + lk * lk - li * li) / (2 * lj * lk)); };
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
        double theta1 = cosine_law(he->edge()->length(), he_next->edge()->length(), he_prev->edge()->length());
        he->edge()->prop("weight") = 0.5 / tan(theta1);
        he->touched() = true;

        CHalfEdge * he_dual = he->he_sym();
        if (!he_dual || he_dual->touched()) continue;
        CHalfEdge * he_dual_next = he_dual->he_next();
        CHalfEdge * he_dual_prev = he_dual->he_prev();
        double theta2 = cosine_law(he_dual->edge()->length(), he_dual_next->edge()->length(), he_dual_prev->edge()->length());
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
    auto fun = [=](double x) { return a * x * x + b * x + c; };
    double f0 = fun(x0);
    double f1 = fun(x1);
    x = x0;
    double fm = f0;
    if(f1 < f0)
    {
        x = x1;
        fm = f1;
    }
    double x2 = -b/a/2.0;
    double f2 = fun(x2);
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
    auto e1 = x->edge;
    auto e2 = y->edge;
    double de = graph->distance(e1, e2, n1, n2);
    return x->length + de + y->length;
}


int CGraphHarmoicMap::calculateBarycenter(CVertex * v)
{
    vector<CVertex*> nei;
    vector<double> ew;
    for (VertexOutHalfedgeIterator heit(mesh, v); !heit.end(); ++heit)
    {
        CHalfEdge * he = *heit;
        nei.push_back(he->target());
        ew.push_back(he->edge()->prop("weight"));
    }
    double a = 0.0;
    for (double w : ew) a += w;
    
    vector<double> fm;
    vector<CTarget*> vx;
    for (SmartGraph::EdgeIt e(*graph); e != INVALID; ++e)
    {
        if (e->source == e->target) // if e is a loop
        {
            // find those nei on this loop

            // find minimum point in this loop
        }
        

        else // if e is not loop
        {

            // find those nei on this edge

            // find minimum point on this edge

        }
        
    }

    return 0;
}

/*
 * for every vertex, move it to its neighbors' weighted center
*/
int CGraphHarmoicMap::harmonicMap()
{
    calculateEdgeLength();
    calculateEdgeWeight();
    int k = 0;
    while (k < 100)
    {
        for (MeshVertexIterator vit(mesh); !vit.end(); ++vit)
        {
            calculateBarycenter(*vit);
        }
    }

    return 0;
}

int CGraphHarmoicMap::writeMap(string filename)
{
    return 0;
}
