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

class mTriangle{
public:	
	ofPoint p[3];
	ofPoint n;
	mTriangle() {		
		for (int i = 0; i < 3; i++) { p[i] = ofPoint(0,0,0); }		
		n.set(0, 1, 0);
	}		
	mTriangle(const mTriangle& t) {
		n = t.n;
		for (int i = 0; i < 3; i++) { p[i] = t.p[i]; }
	}
	mTriangle operator=(const mTriangle& t) {
		n = t.n;
		for (int i = 0; i < 3; i++) { p[i] = t.p[i]; }
	}
		
	void calcnormal(bool invertnormals = false, bool normalize = false) {
		ofPoint nv;
		ofPoint v1, v2;					
		
		v1 = p[1] - p[0];
		v2 = p[2] - p[0];
			
			nv.x = (v1.y * v2.z) - (v1.z * v2.y);
			nv.y = (v1.z * v2.x) - (v1.x * v2.z);
			nv.z = (v1.x * v2.y) - (v1.y * v2.x);

		if (!invertnormals) {
				n.x = nv.x;
				n.y = nv.y;
				n.z = nv.z;		
		} else {
				n.x = -nv.x;
				n.y = -nv.y;
				n.z = -nv.z;
		}

		if(normalize) n.normalize();
		
//		// vertex normals weighted by inverse tri area		
//		float onearea = 1.0f / ((v1.length()*v2.length())*0.5);
//		n *= onearea;
	}
};