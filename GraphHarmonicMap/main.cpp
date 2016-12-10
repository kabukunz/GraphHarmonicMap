// GraphHarmonicMap.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GraphHarmoicMap.h"

int main()
{
    CGraphHarmoicMap * map = new CGraphHarmoicMap();
    map->setMesh("eight.5.m");
    map->setGraph("eight.target");
    //map->test();
    map->harmonicMap();
	map->writeMap("eight.harmonic");

    return 0;
}

