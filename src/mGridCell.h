/*
 
 oF marching cubes implementation by andre sier
 http://s373.net/code/marchingcubes/
 
 Marching cubes implementation after Paul Bourke polygonize voxel.
 http://paulbourke.net/geometry/polygonise/
 
 201207 - of version
 201004 - processing version
 
 code released under http://www.gnu.org/licenses/lgpl-3.0.txt
 copyright 2012, Andre Sier
 http://s373.net
 
 
 */


#pragma once

#include "ofMain.h"

class mGridCell {
public:
	ofPoint p[8];
	float val[8];
	
	mGridCell () {
		for (int i=0; i<8;i++) {
			p[i] = ofPoint(0,0,0); 
			val[i] = 0.0f;
		}
	}
};