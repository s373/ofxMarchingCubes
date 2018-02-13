/*

 oF marching cubes implementation by andre sier
 Marching cubes implementation after Paul Bourke polygonize voxel.
 http://paulbourke.net/geometry/polygonise/

 20161102 014 - adding marching tetrahedron
 201207 - of version
 201004 - processing version

 code released under http://www.gnu.org/licenses/lgpl-3.0.txt
 copyright 2012, s373.net/x art codex studios

 */


#pragma once

#include "ofMain.h"
#include "mTables.h"
#include "mTriangle.h"
#include "mGridCell.h"
#include <vector>
#include <iostream>
#include <fstream>

// using namespace std;

class ofxMarchingCubes;
class ofxMarchingCubes {
public:
	ofxMarchingCubes();
	virtual ~ofxMarchingCubes();

	enum BOUNDMODE{
		ADD = 0,
		CLAMP,
		WRAP
	};

	BOUNDMODE boundmode;
	void setBoundMode(int bm){
		boundmode = (BOUNDMODE) bm;
	}

	inline void addDataIdx(int idx, float val){

		data[idx] += val;

		switch(boundmode){
			case ADD	:  default:  break;
			case CLAMP	:  data[idx]<0 ? data[idx] = 0 : data[idx]>1 ? data[idx]=1.0f: true; break;
			case WRAP	:  data[idx]<0 ? data[idx]+= 1 : data[idx]>1 ? data[idx]-=1: true; break;
		}
	}

	inline void setData(int idx, const float val){
		data[idx]=val;
	}


	void setup(const float sx, const float sy, const float sz,
				const int x, const int y,const  int z);
	void setWorldDim(const float x,const  float y, const float z);
	void setWorldDim(const ofPoint &dim);
	void initResolution(const int x, const int y, const int z);
	string& getinfo();

	void clear();
	void addIsoPoint(const ofPoint &ptpos, const float ptval);
	void addIsoPoint(const int idx, const float ptval);
	void addCube(const ofPoint &pos, const ofPoint &dim, const float ptval);
	void addCube(const int centerx, const int centery, const int centerz,
				 const int dimx, const int dimy, const int dimz, const float ptval);
	void addBall(const ofPoint &pos, const ofPoint &dim, const float ptforce);
	void addBall(const int centerx, const int centery, const int centerz,
				 const int dimx, const int dimy, const int dimz, const float ptval);
	void zeroData();
	void setRndData(const float mn, const float mx);
	void setProbRndData(const float prob,const float mn, const float mx);
	void setRndData(const float mx);


	inline int getDataSize(){
		return data.size();
	}
	inline int getNumTriangles(){
		return trilist.size();
	}

	void setData(const vector <float> &d);
	void addData(const vector <float> &d);
	void multData(const float d);
	void normalizeDataTo(const float val);
	void dataInvert();
	void dataSubs(const float v);
	vector <float> & getData();

	void setDataXY(const vector <float> &d);
	void addDataXY(const vector <float> &d);
	void setDataZSpeed(const float speed);
	void updateDataZspeed();

	void polygoniseData();
	bool isEmpty();
	void checkMinMax();

	virtual void draw(ofPolyRenderMode drawMode=OF_MESH_FILL);
	void toMesh(); // assume load vertlist before
	// void drawwireframe();
	void drawnormals(const float s);

	void calcNormals();
	void setInvertNormals(bool b);
	inline void setOm(int o){om = o;}
	inline void setIso(float  i){isolevel = i;}
	inline float getIso(){return isolevel;}

	virtual void readStl(const string &fn);
	virtual void saveStl(const string &fn);
	virtual void readBinaryStl(const string &fn);
	virtual void saveBinaryStl(const string &fn);
//	void readAsciiStl(string fn);
//	void saveAsciiStl(string fn);

protected:
	inline int Polygonize(const mGridCell & grid, const  float isolevel, vector <mTriangle> &triangles);
	inline int PolygoniseTri(const mGridCell & grid,const float isolevel, vector <mTriangle> &triangles, int v0,int v1,int v2,int v3);
	 inline int getIndex(const ofPoint &ptpos);
	inline ofVec3f VertexInterp(const float isolevel, const ofVec3f &p1, const ofVec3f &p2, const float valp1, const float valp2);

	mGridCell grid;
	float isolevel;
	vector <mTriangle> triangles;
	vector <float> data;
	vector <mTriangle> trilist;
	int ntri;
	int gx, gy, gz, numxyz, gxgy;
	float themin, themax;
	bool invertnormals, normalize;
	bool closesides;
	int om; //0 mc, 1 tetra

	ofPoint worlddim, worldstride, worldcenter, datastride;
	string info, header;

	int	mcZhead, mcZMax;
	float mcZheadf, mcZheadSpeed;


};
