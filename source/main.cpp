// GraphHarmonicMap.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GraphHarmonicMap.h"
#include <time.h>

int GraphHarmonicMap(string meshfilename, string graphfilename, string cutfilename, string outfilename, string options)
{
    CGraphHarmonicMap * map = new CGraphHarmonicMap();

    map->setMesh(meshfilename);
    map->setGraph(graphfilename, cutfilename);
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
    //map->test();

    map->output(outfilename);

    return 0;
}

int Decompose(string meshfilename, string graphfilename, string cutfilename, string outfilename, string options)
{
    CGraphHarmonicMap * map = new CGraphHarmonicMap();

    map->setMesh(meshfilename);
    map->setGraph(graphfilename, cutfilename);
    map->initialMap("continue");

    map->decompose();

    map->output(outfilename);

    return 0;
}

int test()
{
    CMesh * mesh = new CMesh();
    mesh->read_m("alex10m.m");
    clock_t start = clock();
    for (MeshFaceIterator fit(mesh); !fit.end(); ++fit)
    {
        CFace * f = *fit;
        f->prop("area") = 0.3;
        f->prop("normal") = 3;
        f->prop("v1") = 1;
        f->prop("v2") = 2;
        f->prop("v3") = 3;
    }
    clock_t finish = clock();
    cout << "elapsed time is " << (finish - start) << "clicks" << endl;

    double ma = 0;
    for (MeshFaceIterator fit(mesh); !fit.end(); ++fit)
    {
        CFace * f = *fit;
        double fa = f->prop("area");
        if (fa > ma) ma = fa;
    }
    cout << "max area = " << ma << endl;
    return 0;
}

int main(int argc, char * argv[])
{
    string command(argv[1]);
    vector<int> ind;
    string meshfilename;
    string graphfilename;
    string cutfilename;
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
        cout << "usage: GraphHarmonicMap [harmonic|decompose] meshfilename graphfilename [cutfilename] [outfilename] [--options]" << endl;
        return -1;
    }

    meshfilename = string(argv[ind[0]]);
    graphfilename = string(argv[ind[1]]);
    if (ind.size() >= 3)
    {
        cutfilename = string(argv[ind[2]]);
    }
    if (ind.size() >= 4)
    {
        outfilename = string(argv[ind[3]]);
    }

    if (command == "harmonic")
    {
        GraphHarmonicMap(meshfilename, graphfilename, cutfilename, outfilename, options);
    }

    if (command == "decompose")
    {
        Decompose(meshfilename, graphfilename, cutfilename, outfilename, options);
    }

    return 0;
}

