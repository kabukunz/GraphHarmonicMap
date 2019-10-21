#include <cstdlib>
#include <queue>
#include <set>
#include <ctime>
#include <cmath>

#ifndef __clang__
#include <omp.h>
#endif

#include "GraphHarmonicMap.h"
#include "parser/parser.h"
#include "Eigen/Eigen"

const double PI = 3.141592653589793;

CGraphHarmonicMap::CGraphHarmonicMap()
{
    mesh = new CMesh();
    graph = new CGraph();
}


CGraphHarmonicMap::~CGraphHarmonicMap()
{
}

int CGraphHarmonicMap::setMesh(const string & filename)
{
    mesh->read_m(filename);
    return 0;
}

/*
* read a target graph, place source vertex randomly on the graph: arc is chosen randomly, always place target at the middle of each arc
*/
int CGraphHarmonicMap::setGraph(const string & graphfilename)
{
    if (!mesh)
    {
        cerr << "read mesh first" << endl;
        exit(-1);
    }

    ifstream graphfile(graphfilename);
    if (graphfile.good())
    {
        string line;
        while (getline(graphfile, line))
        {
            std::istringstream iss(line);
            string type;
            int cid, x, y;
            double len;
            iss >> type;
            if (type.length() == 0 || type[0] == '#') continue;
            if (type == "Cut")
            {
                iss >> cid >> x >> y >> len;
                auto edge = graph->addEdge(x, y, len);
                vector<CVertex*> cut;
                int vid = -1;
                while (iss >> vid)
                {
                    cut.push_back(mesh->vertex(vid));
                }
                if(cut.empty())
                {
                    cerr << "cut can not be empty" << endl;
                    exit(-1);
                }
                cuts[cid] = cut;
            }
            if (type == "Pants")
            {
                int sid, seed;
                iss >> sid >> seed;
                seeds[sid] = mesh->vertex(seed);
            }
        }
    }
    else
    {
        cout << "can't open graph file: " << graphfilename << endl;
        exit(-1);
    }

    for (CVertex * v : mesh->vertices())
    {
        v->fixed() = false;
        v->critical() = false;
        v->critical2() = false;
        v->cut() = false;
        v->cut2() = false;
        v->touched() = false;
    }

    for (auto c : cuts)
    {
        int id = c.first;
        auto cut = c.second;
        auto e = graph->g.edgeFromId(id);
        auto u = graph->g.u(e);
        auto v = graph->g.v(e);
        bool isfixed = graph->nodeValence[u] == 1 || graph->nodeValence[v] == 1;
        if (isfixed)
        {
            for (auto vi : cut)
            {
                vi->fixed() = true;
            }
        }
    }
    graph->calculateNodeDistance();

    return 0;
}

int CGraphHarmonicMap::calculateEdgeLength()
{
    for (CEdge * e : mesh->edges())
    {
        CVertex * v1 = e->halfedge(0)->source();
        CVertex * v2 = e->halfedge(0)->target();
        e->length() = (v1->point() - v2->point()).norm();
    }
    return 0;
}

int CGraphHarmonicMap::calculateEdgeWeight()
{
    auto inverse_cosine_law = [](double li, double lj, double lk) { return acos((lj * lj + lk * lk - li * li) / (2 * lj * lk)); };
    for (CEdge * e : mesh->edges())
    {
        CHalfEdge * he0 = e->halfedge(0);
        CHalfEdge * he1 = e->halfedge(1);
        if (he0)
        {
            CHalfEdge * he_next = he0->next();
            CHalfEdge * he_prev = he0->prev();
            double theta = inverse_cosine_law(e->length(), he_next->edge()->length(), he_prev->edge()->length());
            e->weight() = 0.5 / tan(theta);
        }
        if (he1)
        {
            CHalfEdge * he_next = he1->next();
            CHalfEdge * he_prev = he1->prev();
            double theta = inverse_cosine_law(e->length(), he_next->edge()->length(), he_prev->edge()->length());
            e->weight() += 0.5 / tan(theta);
        }
        if (e->weight() < 0)
        {
            e->weight() = 0;
        }
    }
    return 0;
}

int CGraphHarmonicMap::runRicciFlow()
{
    int nv = mesh->num_vertices();
    int ne = mesh->num_edges();
    int nf = mesh->num_faces();

    typedef Eigen::ArrayXi Index;
    typedef Eigen::ArrayXXi IndexX;
    typedef Eigen::ArrayXd Array;
    typedef Eigen::ArrayXXd ArrayX;

    auto cosine_law = [&](Array & li, Array & lj, Array & lk) {return acos((lj*lj + lk*lk - li*li) / (2 * lj*lk)); };

    int k = 0;
    for (CVertex * v : mesh->vertices()) v->index() = k++;

    IndexX face = IndexX::Zero(nf, 3);
    k = 0;
    for (CFace * f : mesh->faces())
    {
        f->index() = k++;
        int i = 0;
        for (CVertex * v : f->vertices())
        {
            face(f->index(), i++) = v->index();
        }
    }
    IndexX edge = IndexX::Zero(ne, 2);
    IndexX eif = -IndexX::Ones(ne, 2);
    k = 0;
    for (CEdge * e : mesh->edges())
    {
        int id = k++;
        e->index() = id;
        edge(id, 0) = e->vertex1()->index();
        edge(id, 1) = e->vertex2()->index();
        CFace * f1 = e->face1();
        CFace * f2 = e->face2();
        if (f1) eif(id, 0) = f1->index();
        if (f2) eif(id, 1) = f2->index();
    }

    Array u = Array::Zero(nv);
    Array el = Array::Zero(ne);
    Array ew = Array::Zero(ne);
    Array vkt = Array::Zero(nv);
    Array vk = Array::Zero(nv);

    for (CVertex * v : mesh->vertices())
    {
        int id = v->index();
        bool critical = v->critical();
        bool fixed = v->fixed();
        if (critical)
        {
            cout << "critical point: " << v->index() << endl;
            vkt(id) = -PI;
        }
        else if (fixed && !v->boundary())
        {
            cout << "singular point: " << v->index() << endl;
            vkt(id) = 2 * PI;
        }
    }
    boundary = new CBoundary(mesh);
    vector<CLoop*> & loops = boundary->loops();
    Array bk = 2 * PI * Array::Ones(nv);
    for (auto loop : loops)
    {
        auto hes = loop->halfedges();
        for (CHalfEdge * he : hes)
        {
            CVertex * v = he->target();
            bk[v->index()] = PI;
        }
    }

    k = 0;
    while (k++ < 20)
    {
        // calculate vertex curvature
        for (int i = 0; i < ne; ++i)
        {
            el(i) = exp(u(edge(i, 0))) + exp(u(edge(i, 1)));
        }
        ArrayX r = ArrayX::Zero(nf, 3);
        ArrayX ca = ArrayX::Zero(nf, 3);
        for (int i = 0; i < nf; ++i)
        {
            r(i, 0) = exp(u(face(i, 0)));
            r(i, 1) = exp(u(face(i, 1)));
            r(i, 2) = exp(u(face(i, 2)));
        }
        Array l0 = r.col(1) + r.col(2);
        Array l1 = r.col(2) + r.col(0);
        Array l2 = r.col(0) + r.col(1);
        ca.col(0) = cosine_law(l0, l1, l2);
        ca.col(1) = cosine_law(l1, l2, l0);
        ca.col(2) = cosine_law(l2, l0, l1);
        vk = bk;
        for (int i = 0; i < nf; ++i)
        {
            int f0 = face(i, 0);
            int f1 = face(i, 1);
            int f2 = face(i, 2);
            vk(f0) -= ca(i, 0);
            vk(f1) -= ca(i, 1);
            vk(f2) -= ca(i, 2);
        }

        double err = abs(vk - vkt).maxCoeff();
        cout << "current error is " << err << endl;
        if (err < EPS || std::isnan(err)) break;

        // calculate edge weight
        Array w = sqrt(r.col(0)*r.col(1)*r.col(2) / (r.col(0) + r.col(1) + r.col(2)));
        ew.setZero();
        for (int i = 0; i < ne; ++i)
        {
            if (eif(i, 0) >= 0) ew(i) += w(eif(i, 0)) / el(i);
            if (eif(i, 1) >= 0) ew(i) += w(eif(i, 1)) / el(i);
        }

        // Newton's method
        Eigen::VectorXd b = vkt - vk;
        typedef Eigen::Triplet<double> T;
        vector<T> triplets;
        triplets.reserve(ne * 2 + nv);
        Eigen::SparseMatrix<double> A(nv, nv);
        Eigen::VectorXd x;
        for (int i = 0; i < ne; ++i)
        {
            int v0 = edge(i, 0);
            int v1 = edge(i, 1);
            triplets.push_back(T(v0, v1, -ew(i)));
            triplets.push_back(T(v1, v0, -ew(i)));
            triplets.push_back(T(v0, v0, ew(i)));
            triplets.push_back(T(v1, v1, ew(i)));
        }
        triplets.push_back(T(0, 0, 1));
        A.setFromTriplets(triplets.begin(), triplets.end());
        A.finalize();

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success)
        {
            cerr << "waring: Eigen decomposition failed" << endl;
        }
        x = solver.solve(b);
        Array xa = x.array();
        u += xa - x.mean();
    }
    for (CVertex * v : mesh->vertices())
    {
        v->u() = u[v->index()];
        ostringstream oss;
        oss << " u=(" << v->u() << ")";
        v->string() += oss.str();
    }
    for (CEdge * e : mesh->edges())
    {
        e->length() = el(e->index());
        ostringstream oss;
        oss << " length=(" << e->length() << ")";
        e->string() += oss.str();
    }

    return 0;
}

/*
 * find the minimum point of a quadratic function in interval [x0, x1]
 */
inline double quadraticMininum(double a, double b, double c, double x0, double x1, double &x)
{
    auto func = [=](double x) { return a * x * x + b * x + c; };

    double f0 = func(x0);
    double f1 = func(x1);
    x = x0;
    double fm = f0;
    if (f1 < f0)
    {
        x = x1;
        fm = f1;
    }
    double x2 = -b / a / 2.0;
    double f2 = func(x2);
    if (x0 < x2 && x2 < x1)
    {
        if (f2 < fm)
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
double CGraphHarmonicMap::distance(CTarget * x, CTarget * y)
{
    // first compute the shorest distance between the node of two edges
    auto & ex = x->edge;
    auto & ey = y->edge;
    //double elx = graph->edgeLength[ex];
    double ely = graph->edgeLength[ey];
    //const auto & ux = graph->g.u(ex);
    //const auto & vx = graph->g.v(ex);
    const auto & uy = graph->g.u(ey);
    const auto & vy = graph->g.v(ey);
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
            if (ny == uy) dey = dey + y->length;
            else dey += ely - y->length;
        }

    }
    return dey;
}

double CGraphHarmonicMap::distance(CTarget * x, const SmartGraph::Edge & e, SmartGraph::Node & nx, SmartGraph::Node & ne)
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

double CGraphHarmonicMap::distance(CTarget * x, const SmartGraph::Node & n, SmartGraph::Node & nx)
{
    //bool isLoop = false;
    auto e = x->edge;
    double el = graph->edgeLength[e];

    auto u = graph->g.u(e);
    auto v = graph->g.v(e);



    if (u == v)
    {
        double dx = x->length;
        if (dx > el / 2.0) dx = el - dx;
        double d = graph->distance(n, u);
        nx = u;
        return dx + d;
    }
    else
    {
        double d1 = graph->distance(n, u);
        double d2 = graph->distance(n, v);
        if (x->node == u)
        {
            double dx1 = x->length + d1;
            double dx2 = el - x->length + d2;
            if (dx1 < dx2)
            {
                nx = u;
                return dx1;
            }
            else
            {
                nx = v;
                return dx2;
            }
        }
        else
        {
            double dx1 = el - x->length + d1;
            double dx2 = x->length + d2;
            if (dx1 < dx2)
            {
                nx = u;
                return dx1;
            }
            else
            {
                nx = v;
                return dx2;
            }
        }
    }
}

inline double CGraphHarmonicMap::distance(CTarget * x, SmartGraph::Node n)
{
    SmartGraph::Node nx;
    return distance(x, n, nx);
}


double CGraphHarmonicMap::calculateBarycenter(CVertex * v)
{
    bool isfixed = v->fixed();
    if (isfixed) return 0.0;

    CTarget* * neit = v->neit();
    double * ew = v->ew();
    double a = v->ewsum();
    int nn = v->nn();
    SmartGraph::Node dummyNode;
    double * bx = v->bx();
    short  * be = v->be();
    double * bp = v->bp();
    CTarget* * vx = v->vx();
    double * fm = v->fm();
    
    for (SmartGraph::EdgeIt e(graph->g); e != INVALID; ++e)
    {
        auto ue = graph->g.u(e);
        auto ve = graph->g.v(e);
        double el = graph->edgeLength[e];
        int eid = graph->g.id(e);
        bp[0] = 0;
        bp[1] = el / 2.0;
        bp[2] = el;
        for (int i = 0; i < nn; ++i)
        {
            CTarget * vt = neit[i];
            short b = vt->edge == e;
            if (b)
            {
                if (vt->length < el / 2.0) bp[i+3] = (vt->length + el / 2.0);
                else bp[i+3] = (vt->length - el / 2.0);
            }
            be[i] = b;
        }
        if (ue == ve) // if e is a loop
        {
            sort(bp, bp+nn+3);
            for (int i = 0; i < nn; ++i)
            {
                if (be[i])
                {
                    bx[i] = -neit[i]->length;
                }
                else
                {
                    bx[i] = distance(neit[i], ue);
                }
            }
            // find minimum point in this loop
            CTarget * t = vx[eid];
            t->edge = e;
            t->node = ue;
            t->length = 0.0;
            double ei = 1.0 / EPS;
            for (int i = 0; i < nn + 2; ++i)
            {
                double x0 = bp[i];
                double x1 = bp[i + 1];
                double xm = (x0 + x1) / 2.0;

                double b = 0.0, c = 0.0;
                for (int j = 0; j < nn; ++j)
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
            fm[eid] = ei;
        }

        else // if e is not loop
        {
            for (int i = 0; i < nn; ++i)
            {
                if (be[i])
                {
                    double x = -neit[i]->length;
                    bx[i] = x;
                }
                else
                {
                    double du = distance(neit[i], ue, dummyNode);
                    double dv = distance(neit[i], ve, dummyNode);
                    if (du <= dv) bx[i] = du;
                    else bx[i] = -dv - el;
                }
            }
            double b = 0.0, c = 0.0;
            for (int j = 0; j < nn; ++j)
            {
                b += 2 * ew[j] * bx[j];
                c += ew[j] * bx[j] * bx[j];
            }
            CTarget * t = vx[eid];
            t->edge = e;
            t->node = ue;
            t->length = 0.0;
            double x = 0;
            double mi = quadraticMininum(a, b, c, 0, el, x);
            t->length = x;
            fm[eid] = mi;
        }

    }

    double fmm = fm[0];
    CTarget * vm = vx[0];
    for (int i = 1; i < graph->g.edgeNum(); ++i)
    {
        if (fm[i] < fmm)
        {
            fmm = fm[i];
            vm = vx[i];
        }
    }
    CTarget * vt = v->target();
    double dv = distance(vm, vt);
    vt->edge = vm->edge;
    vt->node = vm->node;
    vt->length = vm->length;

    return dv;
}

/*
 * initial map consists of ordinary harmonic map from pants to a "Y" graph
 */
int CGraphHarmonicMap::initialMap(string method)
{
    calculateEdgeLength();
    calculateEdgeWeight();

    if (method == "continue")
    {
        vector<SmartGraph::Edge> edges;
        for (SmartGraph::EdgeIt e(graph->g); e != INVALID; ++e)
        {
            edges.push_back(e);
        }

        for (CVertex * v : mesh->vertices())
        {
            string s = v->string();
            MeshLib::CParser parser(s);
            vector<MeshLib::CToken*> & tokens = parser.tokens();

            bool hasTarget = false;
            for (auto pt : tokens)
            {
                if (pt->key() == "target")
                {
                    string str = strutil::trim(pt->value(), "()");
                    istringstream iss(str);
                    int edgeId = 0;
                    int nodeId = 0;
                    double length = 0.0;
                    iss >> edgeId >> nodeId >> length;
                    CTarget * t = new CTarget();
                    t->edge = graph->g.edgeFromId(edgeId);
                    t->node = graph->g.u(t->edge);
                    t->length = length;
                    v->target() = t;
                    hasTarget = true;
                }
            }
            if (!hasTarget)
            {
                cout << "vertex target is needed" << endl;
                exit(-1);
            }
        }
    }
    else
    {
        int i = 0;
        for (CVertex * v : mesh->vertices())
        {
            v->index() = i++;
            v->cut() = false;
            v->cut2() = false;
            v->x() = 0.0;
            v->y() = 0.0;
        }
        // label cuts
        for (auto c : cuts)
        {
            auto vs = c.second;
            for (auto v : vs)
            {
                v->cut() = true;
            }
        }
        // trace pants
        traceAllPants();

        for (auto p : pantss)
        {
            int id = p.first;
            auto pants = p.second;
            auto node = graph->g.nodeFromId(id);
            embedPants(node, pants);
        }
    }

    for (CVertex * v : mesh->vertices())
    {
        CTarget * vt = v->target();
        auto edge = vt->edge;
        auto node = vt->node;
        double length = vt->length;
        if (node != graph->g.u(edge))
        {
            vt->node = graph->g.u(edge);
            vt->length = graph->edgeLength[edge] - length;
        }
    }

    return 0;
}

/*
 * for every vertex, move it to its neighbors' weighted center
*/
int CGraphHarmonicMap::harmonicMap()
{
    vector<CVertex*> vv;
    for (CVertex * v : mesh->vertices())
    {
        vv.push_back(v);
        auto vr = v->vertices();
        int nn = vr.size();
        v->nn() = nn;
        v->bx() = new double[nn];
        v->be() = new short[nn];
        v->bp() = new double[nn + 3];
        v->fm() = new double[graph->g.edgeNum()];
        v->vx() = new CTarget*[graph->g.edgeNum()];
        for (int i = 0; i < graph->g.edgeNum(); ++i)
        {
            v->vx()[i] = new CTarget();
        }
        v->neit() = new CTarget*[nn];
        v->ew() = new double[nn];
        v->ewsum() = 0;
        for (int i = 0; i < nn; ++i)
        {
            v->neit()[i] = vr[i]->target();
            v->ew()[i] = mesh->edge(v, vr[i])->weight();
            v->ewsum() += v->ew()[i];
        }
    }

    time_t start = time(NULL);
    int k = 0;
    while (k < 2000)
    {
        double err = 0;
        //random_shuffle(vv.begin(), vv.end());
        #pragma omp parallel for
        for (int i = 0; i < vv.size(); ++i)
        {
            double d = calculateBarycenter(vv[i]);
            if (d > err) err = d;
        }
        if (k % 100 == 0) cout << "#" << k << ": " << err << endl;
        //if (k % 500 == 0) writeMap("harmonic." + to_string(int(k / 500)));
        if (err < EPS) break;
        ++k;
    }
    time_t finish = time(NULL);
    cout << "elapsed time is " << (finish - start) << "s" << endl;

    for (CVertex * v : mesh->vertices())
    {
        CTarget * vt = v->target();
        int eid = -1;
        int nid = -1;
        double length = 0.0;
        double x = 0.0;
        int sign = 0;
        if (vt)
        {
            auto e = vt->edge;            
            double el = graph->edgeLength[e];
            eid = graph->g.id(e);
            nid = graph->g.id(vt->node);
            length = vt->length;
            x = vt->length;
            if (x > el / 2.0) x = (el - x);
        }
        ostringstream oss;
        oss << "uv=(" << x << " 0.43) target=(" << eid << " " << nid << " " << length << ")";
        v->string() = oss.str();
    }

    for (CVertex * v : mesh->vertices())
    {
        double * bx = v->bx();
        delete[] bx;
        short * be = v->be();
        delete[] be;
        double * bp = v->bp();
        delete[] bp;
        double * fm = v->fm();
        delete[] fm;
        
        for (int i = 0; i < graph->g.edgeNum(); ++i)
        {
            delete v->vx()[i];
        }
        delete[] v->vx();

        delete[] v->neit();
        delete[] v->ew();
    }

    return 0;
}

int CGraphHarmonicMap::traceAllPants()
{
    for (CVertex * v : mesh->vertices())
    {
        v->pants() = -1;
    }
    for (auto s : seeds)
    {
        int id = s.first;
        CVertex * seed = s.second;
        vector<CVertex*> pants;
        tracePants(id, seed, pants);
        pantss[id] = pants;
    }
    for (CVertex * v : mesh->vertices())
    {
        if(v->pants() == -1 && !v->cut())
        {
            cout << v->id() << endl;
        }
        assert(v->pants() != -1 || v->cut());
    }
    return 0;
}

int CGraphHarmonicMap::tracePants(int id, CVertex * seed, vector<CVertex*> & pants)
{
    for (CVertex * v : mesh->vertices()) v->touched() = false;

    queue<CVertex*> qe;
    qe.push(seed);
    pants.clear();
    while (!qe.empty())
    {
        CVertex * v = qe.front();
        qe.pop();
        pants.push_back(v);
        v->touched() = true;
        v->pants() = id;
        bool isCut = v->cut();
        if (isCut) continue;
        for (CVertex * vj : v->vertices())
        {
            if (vj->touched()) continue;
            bool isCut = vj->cut();
            if (!isCut) qe.push(vj);
            else pants.push_back(vj);
            vj->touched() = true;
        }
    }
    return 0;
}

/*
 * embed pants by computing a harmonic map from pants to a "Y" shape graph
 */
int CGraphHarmonicMap::embedPants(SmartGraph::Node & node, vector<CVertex*> & pants)
{
    vector<SmartGraph::Edge> edges;
    for (SmartGraph::OutArcIt oa(graph->g, node); oa != INVALID; ++oa)
    {
        SmartGraph::Edge e = oa;
        edges.push_back(e);
    }

    if (edges.size() != 3)
    {
        cerr << "not a valid pants decomposition" << endl;
        exit(-1);
    }
    // if there is a loop
    if (edges[0] == edges[1])
    {
        return embedPants(node, pants, edges[0], edges[1], edges[2]);
    }
    else if (edges[0] == edges[2])
    {
        return embedPants(node, pants, edges[0], edges[2], edges[1]);
    }
    else if (edges[1] == edges[2])
    {
        return embedPants(node, pants, edges[1], edges[2], edges[0]);
    }
    else
    {
        return embedPants(node, pants, edges[0], edges[1], edges[2]);
    }

    return 0;
}

/*
 * e0 may equal to e1, but not to e2
 */
int CGraphHarmonicMap::embedPants(SmartGraph::Node & node, vector<CVertex*> & pants, SmartGraph::Edge & e0, SmartGraph::Edge & e1, SmartGraph::Edge & e2)
{
    // renumber pants vertex
    int k = 0;
    for (auto v : pants)
    {
        v->index() = k++;
    }

    if (e0 == e1)
    {
        int eid = graph->g.id(e0);
        double length = graph->edgeLength[e0];
        auto & cut = cuts[eid];
        vector<CVertex*> vs1, vs2;
        findNeighbors(cut, vs1, vs2);

        // find two side neighbors of cut
        for (auto v : vs1)
        {
            v->x() = -length / 2.0;
            v->y() = 0.0;
            v->cut2() = true;
            v->touched() = true;
        }
        for (auto v : vs2)
        {
            v->x() = 0.0;
            v->y() = -length / 2.0;
            v->cut2() = true;
            v->touched() = true;
        }
        for (auto v : cut)
        {
            v->x() = -length / 2.0;
            v->y() = -length / 2.0;
            v->cut2() = true;
            v->touched() = true;
        }
    }
    else
    {
        auto u0 = graph->g.u(e0);
        auto v0 = graph->g.v(e0);
        if (u0 != node && v0 != node)
        {
            cerr << "edge does not connect to node, graph configuration is not correct" << endl;
            return -1;
        }
        SmartGraph::Node n0 = u0 == node ? v0 : u0;

        int e0id = graph->g.id(e0);
        auto & cut0 = cuts[e0id];
        double l0 = graph->edgeLength[e0];
        bool on0 = graph->nodeValence[n0] == 1;
        for (auto v : cut0)
        {
            v->x() = on0 ? -l0 : -l0 / 2.0;
            v->y() = 0.0;
            v->cut2() = true;
            v->touched() = true;
        }

        auto u1 = graph->g.u(e1);
        auto v1 = graph->g.v(e1);
        if (u1 != node && v1 != node)
        {
            cerr << "edge does not connect to node, graph configuration is not correct" << endl;
            return -1;
        }
        SmartGraph::Node n1 = u1 == node ? v1 : u1;

        int e1id = graph->g.id(e1);
        auto & cut1 = cuts[e1id];
        double l1 = graph->edgeLength[e1];
        bool on1 = graph->nodeValence[n1] == 1;
        for (auto v : cut1)
        {
            v->x() = 0.0;
            v->y() = on1 ? -l1 : -l1 / 2.0;
            v->cut2() = true;
            v->touched() = true;
        }
    }

    auto u2 = graph->g.u(e2);
    auto v2 = graph->g.v(e2);

    if (u2 != node && v2 != node)
    {
        cerr << "edge does not connect to node, graph configuration is not correct" << endl;
        return -1;
    }

    SmartGraph::Node n2 = u2 == node ? v2 : u2;

    int eid = graph->g.id(e2);
    double length = graph->edgeLength[e2];
    auto & cut = cuts[eid];
    bool on2 = graph->nodeValence[n2] == 1;
    for (auto v : cut)
    {
        v->x() = on2 ? length : length / 2.0;
        v->y() = on2 ? length : length / 2.0;
        v->cut2() = true;
        v->touched() = true;
    }

    // compute harmonic map
    typedef Eigen::Triplet<double> T;
    vector<T> triplets, tripletsb;
    int nv = pants.size();
    triplets.reserve(nv * 2);
    Eigen::SparseMatrix<double> A(nv, nv);
    Eigen::SparseMatrix<double> B(nv, nv);
    Eigen::VectorXd x, y;
    Eigen::VectorXd bx, by;
    bx.setZero(nv, 1);
    by.setZero(nv, 1);

    for (auto v : pants)
    {
        int id = v->index();
        bool isCut = v->cut();
        bool isCut2 = v->cut2();
        bool c = isCut || isCut2;
        if (!c)
        {
            for (CVertex * v2 : v->vertices())
            {
                int id2 = v2->index();
                bool isCut = v2->cut();
                bool isCut2 = v2->cut2();
                bool c2 = isCut || isCut2;
                CEdge * e = mesh->edge(v, v2);
                double w = e->weight();
                triplets.push_back(T(id, id, -w));
                if (!c2)
                {
                    triplets.push_back(T(id, id2, w));
                }
                else
                {
                    tripletsb.push_back(T(id, id2, w));
                }
            }
        }
        else
        {
            triplets.push_back(T(id, id, 1));
            bx[id] = v->x();
            by[id] = v->y();
        }
    }


    A.setFromTriplets(triplets.begin(), triplets.end());
    A.finalize();

    B.setFromTriplets(tripletsb.begin(), tripletsb.end());
    B.finalize();

    bx -= B * bx;
    by -= B * by;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Eigen::Success)
    {
        cerr << "waring: Eigen decomposition failed" << endl;
    }
    x = solver.solve(bx);
    y = solver.solve(by);

    for (auto v : pants)
    {
        int id = v->index();
        double xi = x[id];
        double yi = y[id];
        v->x() = xi;
        v->y() = yi;
    }

    double el1 = graph->edgeLength[e1];
    for (auto v : pants)
    {
        double x = v->x();
        double y = v->y();
        CTarget * t = new CTarget();
        t->node = node;
        if (x >= 0 && y >= 0)
        {
            t->edge = e2;
            t->length = x < y ? x : y;
        }
        else if (x < y)
        {
            t->edge = e0;
            t->length = -x;
        }
        else if (y <= x)
        {
            t->edge = e1;
            if (e0 == e1)
            {
                t->length = el1 + y;
            }
            else
            {
                t->length = -y;
            }
        }
        v->target() = t;
    }
    Cut cut3;
    cut3[graph->g.id(e0)] = cuts[graph->g.id(e0)];
    cut3[graph->g.id(e1)] = cuts[graph->g.id(e1)];
    cut3[graph->g.id(e2)] = cuts[graph->g.id(e2)];
    for (auto c : cut3)
    {
        int id = c.first;
        auto cut = c.second;
        auto e = graph->g.edgeFromId(id);
        auto u = graph->g.u(e);
        auto v = graph->g.v(e);
        SmartGraph::Node n = u == node ? v : u;
        for (auto vi : cut)
        {
            bool isfixed = vi->fixed();
            if (isfixed)
            {
                CTarget * vt = vi->target();
                if(!vt)
                {
                    cerr << "vertex has no target:" << vi->id() << endl;
                    exit(-1);
                }
                vt->edge = e;
                vt->node = n;
                vt->length = 0;
                vi->target() = vt;
            }

        }
    }

    return 0;
}

int CGraphHarmonicMap::findNeighbors(vector<CVertex*> & cut, vector<CVertex*> & vs1, vector<CVertex*> & vs2)
{
    vs1.clear();
    vs2.clear();
    for (int i = 0; i < cut.size() - 1; ++i)
    {
        CVertex * v0 = mesh->vertex(cut[i]->id());
        CVertex * v1 = mesh->vertex(cut[i + 1]->id());
        CEdge * e = mesh->edge(v0, v1);
        CHalfEdge * he0 = e->halfedge(0);
        CHalfEdge * he1 = e->halfedge(1);
        CHalfEdge * he = NULL;
        if (he0->source() == v0 && he0->target() == v1)  he = he0;
        if (he1 && he1->source() == v0 && he1->target() == v1) he = he1;
        if (!he)
        {
            cerr << "cut can not be on boundary" << endl;
            return -1;
        }
        vs1.push_back(he->next()->target());
        vs2.push_back(he->dual()->next()->target());
    }
    return 0;
}



void CGraphHarmonicMap::test()
{

}

int CGraphHarmonicMap::traceCriticalTrajectory()
{
    graph->colorize();

    return 0;
}

int CGraphHarmonicMap::decompose()
{
    int k = 0;
    for (CVertex * v : mesh->vertices()) v->index() = k++;
    k = 0;
    for (CFace * f : mesh->faces()) f->index() = k++;

    for (CVertex * v : mesh->vertices())
    {
        bool isfixed = v->fixed();
        CTarget * vt = v->target();
        SmartGraph::Edge e = vt->edge;
        double el = graph->edgeLength[e];
        double l = vt->length;
        if (l > el*0.75) l = el - l;
        if (!isfixed && l <= el*EPS)
        {
            cout << "critical vertex: " << v->index() << endl;
            v->critical2() = true;
        }
    }

    dmesh = new CDynamicMesh(mesh);
    for (CFace * f : mesh->faces())
    {
        if (hasCriticalPoint(f))
        {
            CVertex * v3 = locateCriticalPoint(f);
            v3->critical() = true;
            v3->critical2() = true;
            v3->fixed() = true;
            cout << "critical point inside triangle: " << v3->index() << endl;
        }
    }

    for (CEdge * e : mesh->edges())
    {
        if (hasCriticalPoint(e))
        {
            CVertex *v2 = locateCriticalPoint(e);
            v2->critical2() = true;
        }
    }

    mesh = dmesh;

    // label face
    for (CFace * f : mesh->faces())
    {
        CHalfEdge * he = f->halfedge();
        CVertex * v0 = he->source();
        CVertex * v1 = he->target();
        CVertex * v2 = he->next()->target();
        CVertex * vs;
        // v0,v1,v2 must have same edge id or been critical
        bool c0 = v0->critical2();
        bool c1 = v1->critical2();
        bool c2 = v2->critical2();
        bool ct0 = v0->target();
        bool ct1 = v1->target();
        bool ct2 = v2->target();
        if (c0 && c1 && c2)
        {
            c0 = ct0;
            c1 = ct1;
            c2 = ct2;
        }
        else
        {
            c0 = !c0 && ct0;
            c1 = !c1 && ct1;
            c2 = !c2 && ct2;
        }
        /*if (c0 && c1 && c2)
        {
            cout << "face can't have 3 critical vertices" << endl;
            cout << "face " << f->index() << " is critical: " << hasCriticalPoint(f) << endl;
            exit(-1);
        }*/

        set<int> eids;
        if (c0) {
            bool f0 = v0->fixed();
            if (!f0)
            {
                CTarget * vt0 = v0->target();
                SmartGraph::Edge e0 = vt0->edge;
                eids.insert(graph->g.id(e0));
                vs = v0;
            }
        }
        if (c1) {
            bool f1 = v1->fixed();
            if (!f1)
            {
                CTarget * vt1 = v1->target();
                SmartGraph::Edge e1 = vt1->edge;
                eids.insert(graph->g.id(e1));
                vs = v1;
            }
        }
        if (c2) {
            bool f2 = v2->fixed();
            if (!f2)
            {
                CTarget * vt2 = v2->target();
                SmartGraph::Edge e2 = vt2->edge;
                eids.insert(graph->g.id(e2));
                vs = v2;
            }
        }
        int id = -1;
        if (eids.size() != 1)
        {
            cerr << "face not been splitted: " << f->index() << endl;
            exit(-1);
        }
        else if (eids.size() == 1)
        {
            id = *eids.begin();
        }
        ostringstream oss;
        oss << "e=(" << id << ")";
        oss << " uv=(";

        CTarget * vst = vs->target();

        if (!vst)
        {
            cout << "face " << f->index() << "has 3 critical vertex" << endl;
            exit(-1);
        }
        else
        {
            SmartGraph::Edge es = vst->edge;
            int sign = graph->edgeSign[es];
            double el = graph->edgeLength[es];
            double x[3] = { 0,0,0 };
            bool critical[3];
            for (int i = 0; i < 3; ++i)
            {
                CVertex * v = he->target();
                bool isCritical = v->critical2();
                if (!isCritical)
                {
                    CTarget * vt = v->target();
                    SmartGraph::Edge e = vt->edge;
                    if (e != es)
                    {
                        cout << "e != es" << endl;
                    }
                    x[i] = vt->length;
                    critical[i] = false;
                }
                else
                {
                    critical[i] = true;
                }
                he = he->next();
            }
            if (max(max(x[0], x[1]), x[2]) > el / 2.0)
            {
                if (critical[0]) x[0] = el;
                if (critical[1]) x[1] = el;
                if (critical[2]) x[2] = el;
            }
            if (fabs(max(max(x[0], x[1]), x[2]) - min(min(x[0], x[1]), x[2])) > el / 2.0)
            {
                cout << x[0] << " " << x[1] << " " << x[2] << endl;
            }
            x[0] *= sign;
            x[1] *= sign;
            x[2] *= sign;
            oss << x[0] << " 0.43 " << x[1] << " 0.43 " << x[2] << " 0.43";
        }
        oss << ")";
        f->string() = oss.str();
    }

    // label sharp edge
    for (CEdge * e : mesh->edges())
    {
        CVertex * v1 = e->vertex1();
        CVertex * v2 = e->vertex2();
        bool c1 = v1->critical2();
        bool c2 = v2->critical2();
        if (c1&&c2)
        {
            e->string() = "sharp=(1)";
        }
        else
        {
            //e->string() = "sharp=(0)";
        }
    }

    for (CVertex * v : mesh->vertices())
    {
        if (v->critical2())
        {
            v->string() = "uv=(0 0.43) target=(-1 -1 0)";
        }
    }

    cout << "calculating flat metric" << endl;
    runRicciFlow();

    return 0;
}

bool CGraphHarmonicMap::hasCriticalPoint(CFace * f)
{
    CHalfEdge * he = f->halfedge();

    CEdge * e0 = he->edge();
    CEdge * e1 = he->next()->edge();
    CEdge * e2 = he->prev()->edge();
    return hasCriticalPoint(e0) && hasCriticalPoint(e1) && hasCriticalPoint(e2);
}

bool CGraphHarmonicMap::hasCriticalPoint(CVertex * v1, CVertex * v2)
{
    CEdge * e = mesh->edge(v1, v2);
    if (e) return hasCriticalPoint(e);
    return false;
}

bool CGraphHarmonicMap::hasCriticalPoint(CEdge * e)
{
    CVertex * v0 = e->vertex1();
    CVertex * v1 = e->vertex2();

    bool f0 = v0->fixed();
    bool f1 = v1->fixed();
    if (f0 || f1) return false;

    // if any vertex is critical, then this face has no critical point inside;
    // edge has critical point iff both vertices are critical, that case has been processed
    bool critical0 = v0->critical2();
    bool critical1 = v1->critical2();
    if (critical0 || critical1) return true;

    CTarget * vt0 = v0->target();
    CTarget * vt1 = v1->target();
    SmartGraph::Edge e0 = vt0->edge;
    SmartGraph::Edge e1 = vt1->edge;
    if (e0 != e1) return true;

    // e0 == e1
    auto u = graph->g.u(e0);
    auto v = graph->g.v(e0);
    if (u != v) return false;// not a loop

    // e0 == e1 is a loop
    double el = graph->edgeLength[e0];
    double l0 = vt0->length;
    double l1 = vt1->length;

    if (fabs(l0 - l1) > el / 2.0) return true;
    return false;
}

CVertex * CGraphHarmonicMap::locateCriticalPoint(CFace * f)
{
    CHalfEdge * he = f->halfedge();
    CVertex * v[3];
    v[0] = he->source();
    v[1] = he->target();
    v[2] = he->next()->target();

    CTarget * vt0 = v[0]->target();
    CTarget * vt1 = v[1]->target();
    CTarget * vt2 = v[2]->target();

    SmartGraph::Edge e0 = vt0->edge;
    SmartGraph::Edge e1 = vt1->edge;
    SmartGraph::Edge e2 = vt2->edge;
    double el0 = graph->edgeLength[e0];
    double el1 = graph->edgeLength[e1];
    double el2 = graph->edgeLength[e2];
    double l0 = vt0->length;
    double l1 = vt1->length;
    double l2 = vt2->length;
    if (l0 > el0*0.75) l0 = el0 - l0;
    if (l1 > el1*0.75) l1 = el1 - l1;
    if (l2 > el2*0.75) l2 = el2 - l2;

    double d[3] = { 0,0,0 };
    d[0] = l1*l2 / (l0*l1 + l1*l2 + l2*l0);
    d[1] = l2*l0 / (l0*l1 + l1*l2 + l2*l0);
    d[2] = l0*l1 / (l0*l1 + l1*l2 + l2*l0);;
    CPoint p = v[0]->point()*d[0] + v[1]->point()*d[1] + v[2]->point()*d[2];

    double md = d[0];
    int mi = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (d[i] > md)
        {
            md = d[i];
            mi = i;
        }
    }
    if (1 - md < 1e-3)
    {
        return v[mi];
    }

    CVertex * v3 = dmesh->splitFace(f);
    v3->point() = p;
    return v3;
}

CVertex * CGraphHarmonicMap::locateCriticalPoint(CEdge * e)
{
    CVertex * v0 = e->vertex1();
    CVertex * v1 = e->vertex2();
    CTarget * vt0 = v0->target();
    CTarget * vt1 = v1->target();

    SmartGraph::Edge e0 = vt0->edge;
    SmartGraph::Edge e1 = vt1->edge;

    double el0 = graph->edgeLength[e0];
    double el1 = graph->edgeLength[e1];
    double l0 = vt0->length;
    double l1 = vt1->length;
    if (l0 > el0*0.5) l0 = el0 - l0;
    if (l1 > el1*0.5) l1 = el1 - l1;
    CPoint p = (v0->point()*l1 + v1->point()*l0) / (l0 + l1);
    double d0 = l0 / (l0 + l1);
    double d1 = l1 / (l0 + l1);
    if (d0 < 1e-2) return v0;
    if (d1 < 1e-2) return v1;
    CVertex * v2 = dmesh->splitEdge(e);
    v2->point() = p;
    return v2;
}

int CGraphHarmonicMap::output(string filename)
{
    mesh->write_m(filename);
    return 0;
}

int GraphHarmonicMap(string meshfilename, string graphfilename, string outfilename, string options)
{
    CGraphHarmonicMap * map = new CGraphHarmonicMap();

    map->setMesh(meshfilename);
    map->setGraph(graphfilename);
    if (options == "init")
    {
        map->initialMap("init");
    }
    else if (options == "continue")
    {
        map->initialMap("continue");
    }
    else
    {
        map->initialMap();
    }

    map->harmonicMap();

    map->output(outfilename);

    return 0;
}

int Decompose(string meshfilename, string graphfilename, string outfilename, string options)
{
    CGraphHarmonicMap * map = new CGraphHarmonicMap();

    map->setMesh(meshfilename);
    map->setGraph(graphfilename);
    map->initialMap("continue");

    map->decompose();

    map->output(outfilename);

    return 0;
}

int showUsage()
{
    cout << "usage: GraphHarmonicMap [-harmonic|-decompose] meshfilename graphfilename outfilename [--options]" << endl;
    return 0;
}

int main(int argc, char * argv[])
{
    if (argc < 2)
    {
        return showUsage();
    }
    string command(argv[1]);
    vector<int> ind;
    string meshfilename;
    string graphfilename;
    string outfilename;
    string options;

    for (int i = 2; i < argc; ++i)
    {
        string str(argv[i]);
        if (str.substr(0, 2) == "--")
        {
            options = str.substr(2);
        }
        else
        {
            ind.push_back(i);
        }
    }

    if (ind.size() < 3)
    {
        return showUsage();
    }

    meshfilename = string(argv[ind[0]]);
    graphfilename = string(argv[ind[1]]);

    if (ind.size() >= 3)
    {
        outfilename = string(argv[ind[2]]);
    }

    if (command == "-harmonic")
    {
        GraphHarmonicMap(meshfilename, graphfilename, outfilename, options);
    }

    if (command == "-decompose")
    {
        Decompose(meshfilename, graphfilename, outfilename, options);
    }

    return 0;
}
