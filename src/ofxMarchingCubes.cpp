

#include "ofxMarchingCubes.h"












ofxMarchingCubes::ofxMarchingCubes(){
	isolevel = 0.0f;
	triangles.clear();
	data.clear();
	trilist.clear();
	ntri=0;
	gx = gy = gz = numxyz = gxgy = 0;
	themin = themax = 0;
	invertnormals = true;
	closesides = true;		
	mcZhead = mcZMax = 0;
	mcZheadf = 0;
	mcZheadSpeed = 0.001;
	header = "s373.net";
	for(int i=0; i<72; i++){
		header += (char)ofRandom(48,90);
	}
}










ofxMarchingCubes::~ofxMarchingCubes(){
	triangles.clear();
	data.clear();
	trilist.clear();
}













void ofxMarchingCubes::setup(const float sx, const float sy, const float sz, 
							 const int x, const int y, const int z){

	initResolution(x, y, z);
	setWorldDim(sx, sy, sz);

	int ntris = 5;
	for (int i = 0; i < ntris; i++) {
		mTriangle tri;
		triangles.push_back( tri );
	}
	isolevel = 0.0025f;

	polygoniseData();
	cout << getinfo() << endl;
}








void ofxMarchingCubes::setWorldDim(const float x,const float y,const  float z){
	worlddim.set(x, y, z);
	worldstride.set(gx / x, gy / y, gz / z);
	datastride.set(x / gx, y / gy, z / gz);
	worldcenter.set(x/2.0f, y/2.0f, z/2.0f);
	
}








void ofxMarchingCubes::setWorldDim(const ofPoint &dim) {
	setWorldDim(dim.x, dim.y, dim.z);
}








void ofxMarchingCubes::initResolution(const int x,const  int y,const  int z){
	gx = x;
	gy = y;
	gz = z;
	gxgy = x * y;
	numxyz = x * y * z;	

	data.clear();
	data.reserve(numxyz);
	for (int i = 0; i < numxyz; i++) {
		data.push_back( 0.0f );
	}
	trilist.clear();
	
	mcZhead = 0;
	mcZMax = gz;
	mcZheadf = 0;
	
}











string& ofxMarchingCubes::getinfo(){
	info = "tris: " + ofToString(ntri) + " volume: " + ofToString(gx) + " " + ofToString(gy)
	+ " " + ofToString(gz) + " cells: " + ofToString(numxyz) + " iso: " + ofToString(isolevel);
	return info;	
}











void ofxMarchingCubes::clear(){
	zeroData();	
}








void ofxMarchingCubes::addIsoPoint(const ofPoint &ptpos, const float ptval){
	int idx = getIndex(ptpos);
	data[idx] += ptval;
}








void ofxMarchingCubes::addIsoPoint(const int idx, const float ptval){
	data[idx] += ptval;
}








void ofxMarchingCubes::addCube(const ofPoint &ptpos, const ofPoint &ptdim, const float ptval){
	// replace this getindex	
	int cx = (int) (ptpos.x * worldstride.x);
	int cy = (int) (ptpos.y * worldstride.y);
	int cz = (int) (ptpos.z * worldstride.z);
	
	if (closesides) {
		if (cx < 1) {
			cx = 1;
		}
		if (cy < 1) {
			cy = 1;
		}
		if (cz < 1) {
			cz = 1;
		}
		if (cx >= gx - 1) {
			cx = gx - 2;
		}
		if (cy >= gy - 1) {
			cy = gy - 2;
		}
		if (cz >= gz - 1) {
			cz = gz - 2;
		}
	} else {
		if (cx < 0) {
			cx = 0;
		}
		if (cy < 0) {
			cy = 0;
		}
		if (cz < 0) {
			cz = 0;
		}
		if (cx >= gx) {
			cx = gx - 1;
		}
		if (cy >= gy) {
			cy = gy - 1;
		}
		if (cz >= gz) {
			cz = gz - 1;
		}
	}
	
	int dimx = (int) (ptdim.x * worldstride.x); 
	int dimy = (int) (ptdim.y * worldstride.y);
	int dimz = (int) (ptdim.z * worldstride.z);
	
	if (dimx < 2) {
		dimx = 2;
	}
	if (dimy < 2) {
		dimy = 2;
	}
	if (dimz < 2) {
		dimz = 2;
	}
	addCube( cx, cy, cz, dimx, dimy, dimz, ptval);
}








void ofxMarchingCubes::addCube( const int centerx, const int centery,const int centerz,
							   const int dimx, const int dimy, int dimz, const float val) {
	int hx = dimx / 2;
	int hy = dimy / 2;
	int hz = dimz / 2;
	
	int sx = centerx - hx;
	int sy = centery - hy;
	int sz = centerz - hz;
	
	if (sx < 1) {
		sx = 1;
	}
	if (sy < 1) {
		sy = 1;
	}
	if (sz < 1) {
		sz = 1;
	}
	
	int tx = MIN(centerx + hx, gx - 1);
	int ty = MIN(centery + hy, gy - 1);
	int tz = MIN(centerz + hz, gz - 1);
	
	int idx = 0;
	for (int i = sx; i < tx; i++) {
		for (int j = sy; j < ty; j++) {
			for (int k = sz; k < tz; k++) {
				idx = i + j * gx + k * gxgy;
				data[idx] += val;
			}
		}
	}
}









void ofxMarchingCubes::addBall(const ofPoint &ptpos, const ofPoint &ptdim, const float ptval){
	// replace this getindex	
	int cx = (int) (ptpos.x * worldstride.x);
	int cy = (int) (ptpos.y * worldstride.y);
	int cz = (int) (ptpos.z * worldstride.z);
	
	if (closesides) {
		if (cx < 1) {
			cx = 1;
		}
		if (cy < 1) {
			cy = 1;
		}
		if (cz < 1) {
			cz = 1;
		}
		if (cx >= gx - 1) {
			cx = gx - 2;
		}
		if (cy >= gy - 1) {
			cy = gy - 2;
		}
		if (cz >= gz - 1) {
			cz = gz - 2;
		}
	} else {
		if (cx < 0) {
			cx = 0;
		}
		if (cy < 0) {
			cy = 0;
		}
		if (cz < 0) {
			cz = 0;
		}
		if (cx >= gx) {
			cx = gx - 1;
		}
		if (cy >= gy) {
			cy = gy - 1;
		}
		if (cz >= gz) {
			cz = gz - 1;
		}
	}
	
	int dimx = (int) (ptdim.x * worldstride.x); 
	int dimy = (int) (ptdim.y * worldstride.y);
	int dimz = (int) (ptdim.z * worldstride.z);
	
	if (dimx < 2) {
		dimx = 2;
	}
	if (dimy < 2) {
		dimy = 2;
	}
	if (dimz < 2) {
		dimz = 2;
	}
	addBall( cx, cy, cz, dimx, dimy, dimz, ptval);
}








void ofxMarchingCubes::addBall( const int centerx, const int centery,const int centerz,
							   const int dimx, const int dimy, int dimz, const float val) {
	int hx = dimx / 2;
	int hy = dimy / 2;
	int hz = dimz / 2;
	
	int sx = centerx - hx;
	int sy = centery - hy;
	int sz = centerz - hz;
	
	if (sx < 1) {
		sx = 1;
	}
	if (sy < 1) {
		sy = 1;
	}
	if (sz < 1) {
		sz = 1;
	}
	
	int tx = MIN(centerx + hx, gx - 1);
	int ty = MIN(centery + hy, gy - 1);
	int tz = MIN(centerz + hz, gz - 1);
	
	int idx = 0;
	ofPoint delta, target, src(centerx*worlddim.x, centery*worlddim.y, centerz*worlddim.z);
	for (int i = sx; i < tx; i++) {
		for (int j = sy; j < ty; j++) {
			for (int k = sz; k < tz; k++) {
				idx = i + j * gx + k * gxgy;
				target.set(i*worlddim.x, j*worlddim.y, k*worlddim.z);
				delta = target-src;
				float fval = val / delta.length();
				data[idx] += fval;
			}
		}
	}
}







//
//
//void ofxMarchingCubes::addCube( const int centerx, const int centery,const int centerz,
//							   const int dimx, const int dimy, int dimz, const float val) {
//	int hx = dimx / 2;
//	int hy = dimy / 2;
//	int hz = dimz / 2;
//	
//	int sx = centerx - hx;
//	int sy = centery - hy;
//	int sz = centerz - hz;
//	
//	if (sx < 1) {
//		sx = 1;
//	}
//	if (sy < 1) {
//		sy = 1;
//	}
//	if (sz < 1) {
//		sz = 1;
//	}
//	
//	int tx = MIN(centerx + hx, gx - 1);
//	int ty = MIN(centery + hy, gy - 1);
//	int tz = MIN(centerz + hz, gz - 1);
//	
//	int idx = 0;
//	for (int i = sx; i < tx; i++) {
//		for (int j = sy; j < ty; j++) {
//			for (int k = sz; k < tz; k++) {
//				idx = i + j * gx + k * gxgy;
//				data[idx] += val;
//			}
//		}
//	}
//}
//
//
//






void ofxMarchingCubes::zeroData() {
	for (int i = 0; i < numxyz; i++) {
		data[i] = 0.0f;
	}	
}








void ofxMarchingCubes::setRndData(const float mx) {
	for (int i = 0; i < numxyz; i++) {
		data[i] = ofRandom(mx);
	}
}









void ofxMarchingCubes::setRndData(const float mn,const float mx) {
	for (int i = 0; i < numxyz; i++) {
		data[i] = ofRandom(mn, mx);
	}
}









void ofxMarchingCubes::setData(const vector <float> &d) {
	for (int i = 0; i < numxyz; i++) {
		data[i] = d[i];
	}
	
}








void ofxMarchingCubes::addData(const vector <float> &d) {
	for (int i = 0; i < numxyz; i++) {
		data[i] += d[i];
	}
	
}








void ofxMarchingCubes::multData(const float d) {
	for (int i = 0; i < numxyz; i++) {
		data[i] *= d;
	}
	
}








void ofxMarchingCubes::normalizeDataTo(const float val) {
	checkMinMax();
	float dist = themax - themin;
	float dst = val / dist;// 1.0f / v;
	multData(dst);
	
}








void ofxMarchingCubes::dataInvert() {
	float max = -1;
	for (int i = 0; i < numxyz; i++) {
		if(data[i]>max){
			max = data[i];
		}
	}
	
	for (int i = 0; i < numxyz; i++) {
		data[i] = max - data[i];
	}
	
}








void ofxMarchingCubes::dataSubs(const float v) {
	for (int i = 0; i < numxyz; i++) {
		data[i] = v - data[i];
	}
	
}








vector <float> & ofxMarchingCubes::getData() {
	return data;
}








void ofxMarchingCubes::setDataXY(const vector <float> &d) {
	
	int locz = mcZhead * gxgy;
	for(int i=0; i<gxgy;i++){
		data[locz+i] = d[i];
	}	
	
}








void ofxMarchingCubes::addDataXY(const vector <float> &d) {
	
	int locz = mcZhead * gxgy;
	for(int i=0; i<gxgy;i++){
		data[locz+i] += d[i];
	}	
	
}








void ofxMarchingCubes::setDataZSpeed(const float speed){
	mcZheadSpeed = speed;	
}








void ofxMarchingCubes::updateDataZspeed(){
	mcZheadf += mcZheadSpeed;
	if(mcZheadf>mcZMax){
		mcZheadf-=mcZMax;
	}
	if(mcZheadf<0){
		mcZheadf+=mcZMax;
	}
	
	mcZhead = mcZheadf;
}








void ofxMarchingCubes::polygoniseData() {
	// Polygonise the grid
	ntri = 0;
	trilist.clear();
	for (int i = 0; i < gx - 1; i++) {
		for (int j = 0; j < gy - 1; j++) {
			for (int k = 0; k < gz - 1; k++) {
				grid.p[0].x = i * datastride.x;
				grid.p[0].y = j * datastride.y;
				grid.p[0].z = k * datastride.z;
				grid.val[0] = data[i + j * gx + k * gxgy];
				grid.p[1].x = (i + 1) * datastride.x;
				grid.p[1].y = j * datastride.y;
				grid.p[1].z = k * datastride.z;
				grid.val[1] = data[i + 1 + j * gx + k * gxgy];
				grid.p[2].x = (i + 1) * datastride.x;
				grid.p[2].y = (j + 1) * datastride.y;
				grid.p[2].z = k * datastride.z;
				grid.val[2] = data[i + 1 + (j + 1) * gx + k * gxgy];
				grid.p[3].x = i * datastride.x;
				grid.p[3].y = (j + 1) * datastride.y;
				grid.p[3].z = k * datastride.z;
				grid.val[3] = data[i + (j + 1) * gx + k * gxgy];
				grid.p[4].x = i * datastride.x;
				grid.p[4].y = j * datastride.y;
				grid.p[4].z = (k + 1) * datastride.z;
				grid.val[4] = data[i + j * gx + (k + 1) * gxgy];
				grid.p[5].x = (i + 1) * datastride.x;
				grid.p[5].y = j * datastride.y;
				grid.p[5].z = (k + 1) * datastride.z;
				grid.val[5] = data[i + 1 + j * gx + (k + 1) * gxgy];
				grid.p[6].x = (i + 1) * datastride.x;
				grid.p[6].y = (j + 1) * datastride.y;
				grid.p[6].z = (k + 1) * datastride.z;
				grid.val[6] = data[i + 1 + (j + 1) * gx + (k + 1) * gxgy];
				grid.p[7].x = i * datastride.x;
				grid.p[7].y = (j + 1) * datastride.y;
				grid.p[7].z = (k + 1) * datastride.z;
				grid.val[7] = data[i + (j + 1) * gx + (k + 1) * gxgy];
				
				int n = Polygonize(grid, isolevel, triangles);
				
				// calc tri norms
				for (int a0 = 0; a0 < n; a0++) {
					triangles[a0].calcnormal(invertnormals);
				}
				
				
				for (int l = 0; l < n; l++) {
					mTriangle tri( triangles[l] );
					trilist.push_back(tri);
				}
				ntri += n;
			}
		}
	}
	// println("Total of triangles: "+ntri);
	
	// Now do something with the triangles ....
	// Here I just write them to a geom file
	// http://local.wasp.uwa.edu.au/~pbourke/dataformats/geom/
	// fprintf(stderr,"Writing triangles ...\n");
	// if ((fptr = fopen("output.geom","w")) == NULL) {
	// fprintf(stderr,"Failed to open output file\n");
	// exit(-1);
	// }
	// for (i=0;i<ntri;i++) {
	// fprintf(fptr,"f3 ");
	// for (k=0;k<3;k++) {
	// fprintf(fptr,"%g %g %g ",tri[i].p[k].x,tri[i].p[k].y,tri[i].p[k].z);
	// }
	// fprintf(fptr,"0.5 0.5 0.5\n"); // colour
	// }
	// fclose(fptr);
	
}








bool ofxMarchingCubes::isEmpty() {
	bool empty = true;
	for(int i=0; i< numxyz; i++){
		if(data[i]>0){
			empty = false;
			break;
		}
	}
	return empty;
	
}








void ofxMarchingCubes::checkMinMax() {
	themin = 1e10;
	themax = -1e10;
	for (int i = 0; i < numxyz; i++) {
		if (data[i] > themax) {
			themax = data[i];
		}
		if (data[i] < themin) {
			themin = data[i];
		}
	}	
	
}





//GLfloat tri[12];


void ofxMarchingCubes::draw() {
	
	vector<float> tripts, normpts;
	for(int i=0; i<trilist.size(); i++){
		normpts.push_back(trilist[i].n.x);
		normpts.push_back(trilist[i].n.y);
		normpts.push_back(trilist[i].n.z);		normpts.push_back(trilist[i].n.x);
		normpts.push_back(trilist[i].n.y);
		normpts.push_back(trilist[i].n.z);		normpts.push_back(trilist[i].n.x);
		normpts.push_back(trilist[i].n.y);
		normpts.push_back(trilist[i].n.z);
		tripts.push_back(trilist[i].p[0].x);
		tripts.push_back(trilist[i].p[0].y);
		tripts.push_back(trilist[i].p[0].z);
		tripts.push_back(trilist[i].p[1].x);
		tripts.push_back(trilist[i].p[1].y);
		tripts.push_back(trilist[i].p[1].z);
		tripts.push_back(trilist[i].p[2].x);
		tripts.push_back(trilist[i].p[2].y);
		tripts.push_back(trilist[i].p[2].z);
	}
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &tripts[0]);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, &normpts[0]);
	glVertexPointer(3, GL_FLOAT, 0, &tripts[0]);
	glDrawArrays(GL_TRIANGLES, 0, (int)trilist.size()*3);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	
	
//	for(int i=0; i<trilist.size(); i++){
//		mTriangle &tri = trilist[i];
//		glBegin(GL_TRIANGLES);
//		glNormal3f( tri.n.x, tri.n.y, tri.n.z );
//		glVertex3f( tri.p[0].x, tri.p[0].y, tri.p[0].z );
//		glVertex3f( tri.p[1].x, tri.p[1].y, tri.p[1].z );
//		glVertex3f( tri.p[2].x, tri.p[2].y, tri.p[2].z );
//		glEnd();
//	}
}







void ofxMarchingCubes::drawnormals(const float s) {
	
	
	
	const float one3 = 1.0f/3.0f;
	
	vector<float> linepts;
	for(int i=0; i<trilist.size(); i++){
		float x = (trilist[i].p[0].x + trilist[i].p[1].x + trilist[i].p[2].x) *one3;
		float y = (trilist[i].p[0].y + trilist[i].p[1].y + trilist[i].p[2].y) *one3;
		float z = (trilist[i].p[0].z + trilist[i].p[1].z + trilist[i].p[2].z) *one3;		
		linepts.push_back(x);
		linepts.push_back(y);
		linepts.push_back(z);
		linepts.push_back(x+s*trilist[i].n.x);
		linepts.push_back(y+s*trilist[i].n.y);
		linepts.push_back(z+s*trilist[i].n.z);
	}
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &linepts[0]);
	glDrawArrays(GL_LINES, 0, (int)trilist.size()*2);
	glDisableClientState(GL_VERTEX_ARRAY);
	
//	for (int i = 0; i < trilist.size(); i++) {
//		mTriangle &tri = trilist[i];
//		float x = (tri.p[0].x + tri.p[1].x + tri.p[2].x) *one3;
//		float y = (tri.p[0].y + tri.p[1].y + tri.p[2].y) *one3;
//		float z = (tri.p[0].z + tri.p[1].z + tri.p[2].z) *one3;
//				
//		glBegin(GL_LINES);
//		glVertex3f( x,y,z );
//		glVertex3f( x + tri.n.x * s, y + tri.n.y * s, z + tri.n.z * s );
//		glEnd();
//		
//	}
	
}






void ofxMarchingCubes::readStl(const string &fn) {
	readBinaryStl(fn);	
}


void ofxMarchingCubes::saveStl(const string &fn) {
	saveBinaryStl(fn);	
}



void ofxMarchingCubes::saveBinaryStl(const string &fn) {
	//		UINT8[80] – Header
	//		UINT32 – Number of triangles
	//		
	//		foreach triangle
	//		REAL32[3] – Normal vector
	//		REAL32[3] – Vertex 1
	//		REAL32[3] – Vertex 2
	//		REAL32[3] – Vertex 3
	//		UINT16 – Attribute byte count
	//		end
	
	ofstream myfile;
	myfile.open (ofToDataPath(fn).c_str(), ios::binary );
	
	myfile.precision(4);
	
	myfile.write(header.c_str(), sizeof(char)*80);
	
	short atribute = 0;
	
	int numt=trilist.size();
	myfile.write((char *)&numt,sizeof(int));
	
	
	for(int i=0; i<trilist.size();i++){		
		mTriangle &t = trilist[i];		
		myfile.write((char *)&t.n.x,sizeof(float));
		myfile.write((char *)&t.n.y,sizeof(float));
		myfile.write((char *)&t.n.z,sizeof(float));		
		myfile.write((char *)&t.p[2].x,sizeof(float));
		myfile.write((char *)&t.p[2].y,sizeof(float));
		myfile.write((char *)&t.p[2].z,sizeof(float));		
		myfile.write((char *)&t.p[1].x,sizeof(float));
		myfile.write((char *)&t.p[1].y,sizeof(float));
		myfile.write((char *)&t.p[1].z,sizeof(float));		
		myfile.write((char *)&t.p[0].x,sizeof(float));
		myfile.write((char *)&t.p[0].y,sizeof(float));
		myfile.write((char *)&t.p[0].z,sizeof(float));		
		myfile.write((char *)&atribute,sizeof(short));		
	}
	
//	myfile.write(header.c_str(), sizeof(char)*80);		
	myfile.close();
	
	cout << "ofxMarchingCubes saved binary " << fn << " tris: "<< trilist.size() <<  endl;
	
	
}

void ofxMarchingCubes::readBinaryStl(const string &fn) {

	string file = ofToDataPath(fn, true);	
	ifstream stl_file(file.c_str(), ios::in | ios::binary);		
	
	if(!stl_file.is_open()){
		cout << "ofxMarchingCubes: could not read file " << file << endl;		
	}else{		
		int	through = 0;		
		char h[80];		
		for(int i=0;i<80;i++){ stl_file.read((char *)&h[i], sizeof(char)); }		
		cout << h << endl;
		
		int numtris;
		stl_file.read((char *)&numtris, sizeof(int));
		
		cout << "ofxMarchingCubes reading binary stl " << fn << " numtris = " << ofToString(numtris) << endl;
		
		trilist.clear();
		
		while ((through < numtris))// && (!stl_file.eof()))
		{
			float	floats[12];
			short	attribute;
			
			for (int i=0; i<12; i++) {
				stl_file.read((char *)&floats[i], sizeof(float));
			}
			stl_file.read((char *)&attribute, sizeof(short));
						
			// add triangles
			mTriangle tri;
			tri.n.x = floats[0];
			tri.n.y = floats[1];
			tri.n.z = floats[2];
			tri.p[2].x = floats[3];
			tri.p[2].y = floats[4];
			tri.p[2].z = floats[5];
			tri.p[1].x = floats[6];
			tri.p[1].y = floats[7];
			tri.p[1].z = floats[8];
			tri.p[0].x = floats[9];
			tri.p[0].y = floats[10];
			tri.p[0].z = floats[11];
			trilist.push_back(tri);
			
			through++;
			
//			if(through%1000==0){
//				cout << "adding tri " << through << " " << floats << ": ";
//				for(int i=0; i<12; i++){
//					cout << " " << floats[i];
//				}
//				cout << endl;
//			}
			
		}
		
		stl_file.close();
		
		
		
	}
	
	

	
}




// void ofxMarchingCubes::readAsciiStl(string fn) {
// 	// phps binary is enough
	
// 	ifstream stl_file (fn.c_str(), ios::in );

	

// 	trilist.clear();
	
// 	while (!stl_file.eof()) {
// 		float floats[12];
		
// 	}
	
// //	for(int i=0; i<trilist.size();i++){
// //		
// //		mTriangle &t = trilist[i];
// //		myfile << "facet normal "<< t.n.x << " " << t.n.y << " " << t.n.z << "\n";
// //		myfile << "outer loop" << "\n";
// //		myfile << "vertex " << t.p[0].x << " " << t.p[0].y << " " << t.p[0].z << "\n";
// //		myfile << "vertex " << t.p[1].x << " " << t.p[1].y << " " << t.p[1].z << "\n";
// //		myfile << "vertex " << t.p[2].x << " " << t.p[2].y << " " << t.p[2].z << "\n";
// //		myfile << "end loop" << "\n";
// //		myfile << "endfacet" << "\n";
// //	}
// //	
// //	myfile << "endsolid "<<name << "\n";
// //	myfile.close();
// //	
// //	cout << "ofxMarchingCubes saved ascii " << fn << "tris: "<< trilist.size() <<  endl;
	
	
	
// }





//void ofxMarchingCubes::saveAsciiStl(string fn) {
//
//
//	ofstream myfile;
//	myfile.open (ofToDataPath(fn).c_str());
//	
//	string name = "s373-" + fn;
//	
//	myfile << "solid "<<name << "\n";
//
//	for(int i=0; i<trilist.size();i++){
//	
//		mTriangle &t = trilist[i];
//		myfile << "facet normal "<< t.n.x << " " << t.n.y << " " << t.n.z << "\n";
//		myfile << "outer loop" << "\n";
//		myfile << "vertex " << t.p[0].x << " " << t.p[0].y << " " << t.p[0].z << "\n";
//		myfile << "vertex " << t.p[1].x << " " << t.p[1].y << " " << t.p[1].z << "\n";
//		myfile << "vertex " << t.p[2].x << " " << t.p[2].y << " " << t.p[2].z << "\n";
//		myfile << "end loop" << "\n";
//		myfile << "endfacet" << "\n";
//	}
//	
//	myfile << "endsolid "<<name << "\n";
//	myfile.close();
//	
//	cout << "ofxMarchingCubes saved ascii " << fn << "tris: "<< trilist.size() <<  endl;
//	
//
//	
//}







	
 int	ofxMarchingCubes::Polygonize(const mGridCell grid,const float isolevel,  vector <mTriangle> &triangles) {
	int i, ntriang;
	int cubeindex;
	ofPoint vertlist[12];// = new ofPoint[12];
	/*
	 * Determine the index into the edge table which tells us which vertices
	 * are inside of the surface
	 */
	cubeindex = 0;
	if (grid.val[0] < isolevel) {
		cubeindex |= 1;
	}
	if (grid.val[1] < isolevel) {
		cubeindex |= 2;
	}
	if (grid.val[2] < isolevel) {
		cubeindex |= 4;
	}
	if (grid.val[3] < isolevel) {
		cubeindex |= 8;
	}
	if (grid.val[4] < isolevel) {
		cubeindex |= 16;
	}
	if (grid.val[5] < isolevel) {
		cubeindex |= 32;
	}
	if (grid.val[6] < isolevel) {
		cubeindex |= 64;
	}
	if (grid.val[7] < isolevel) {
		cubeindex |= 128;
	}
	
	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0) {
		return (0);
	}
	
	/* Find the vertices where the surface intersects the cube */
	// int temp = edgeTable[cubeindex] & 1;
	if ((edgeTable[cubeindex] & 1) != 0) {
		vertlist[0] = VertexInterp(isolevel, grid.p[0], grid.p[1],
								   grid.val[0], grid.val[1]);
	}
	if ((edgeTable[cubeindex] & 2) != 0) {
		vertlist[1] = VertexInterp(isolevel, grid.p[1], grid.p[2],
								   grid.val[1], grid.val[2]);
	}
	if ((edgeTable[cubeindex] & 4) != 0) {
		vertlist[2] = VertexInterp(isolevel, grid.p[2], grid.p[3],
								   grid.val[2], grid.val[3]);
	}
	if ((edgeTable[cubeindex] & 8) != 0) {
		vertlist[3] = VertexInterp(isolevel, grid.p[3], grid.p[0],
								   grid.val[3], grid.val[0]);
	}
	if ((edgeTable[cubeindex] & 16) != 0) {
		vertlist[4] = VertexInterp(isolevel, grid.p[4], grid.p[5],
								   grid.val[4], grid.val[5]);
	}
	if ((edgeTable[cubeindex] & 32) != 0) {
		vertlist[5] = VertexInterp(isolevel, grid.p[5], grid.p[6],
								   grid.val[5], grid.val[6]);
	}
	if ((edgeTable[cubeindex] & 64) != 0) {
		vertlist[6] = VertexInterp(isolevel, grid.p[6], grid.p[7],
								   grid.val[6], grid.val[7]);
	}
	if ((edgeTable[cubeindex] & 128) != 0) {
		vertlist[7] = VertexInterp(isolevel, grid.p[7], grid.p[4],
								   grid.val[7], grid.val[4]);
	}
	if ((edgeTable[cubeindex] & 256) != 0) {
		vertlist[8] = VertexInterp(isolevel, grid.p[0], grid.p[4],
								   grid.val[0], grid.val[4]);
	}
	if ((edgeTable[cubeindex] & 512) != 0) {
		vertlist[9] = VertexInterp(isolevel, grid.p[1], grid.p[5],
								   grid.val[1], grid.val[5]);
	}
	if ((edgeTable[cubeindex] & 1024) != 0) {
		vertlist[10] = VertexInterp(isolevel, grid.p[2], grid.p[6],
									grid.val[2], grid.val[6]);
	}
	if ((edgeTable[cubeindex] & 2048) != 0) {
		vertlist[11] = VertexInterp(isolevel, grid.p[3], grid.p[7],
									grid.val[3], grid.val[7]);
	}
	
	/* Create the triangle */
	ntriang = 0;
	for (i = 0; triTable[cubeindex][i] != -1; i += 3) {
		triangles[ntriang].p[0] = vertlist[triTable[cubeindex][i]];
		triangles[ntriang].p[1] = vertlist[triTable[cubeindex][i + 1]];
		triangles[ntriang].p[2] = vertlist[triTable[cubeindex][i + 2]];
		ntriang++;
	}
	
	return (ntriang);
	
}


/*
 * Linearly interpolate the position where an isosurface cuts an edge
 * between two vertices, each with their own scalar value
 */
ofVec3f ofxMarchingCubes::VertexInterp(const float isolevel,const  ofVec3f &p1,const  ofVec3f &p2,const  float valp1,const float valp2){

	float mu;
	ofPoint p(0,0,0);
	float ep = 0.00001;
	
	if (ABS(isolevel - valp1) < ep) {
		return (p1);
	}
	if (ABS(isolevel - valp2) < ep) {
		return (p2);
	}
	if (ABS(valp1 - valp2) < ep) {
		return (p1);
	}
	mu = ((isolevel - valp1) / (valp2 - valp1));
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);
	
	return (p);
}






 int ofxMarchingCubes::getIndex(const ofPoint &ptpos) {
	int cx = (int) (ptpos.x * worldstride.x); 
	int cy = (int) (ptpos.y * worldstride.y);
	int cz = (int) (ptpos.z * worldstride.z);
	if (closesides) {
		if (cx < 1) {
			cx = 1;
		}
		if (cy < 1) {
			cy = 1;
		}
		if (cz < 1) {
			cz = 1;
		}
		if (cx >= gx - 1) {
			cx = gx - 2;
		}
		if (cy >= gy - 1) {
			cy = gy - 2;
		}
		if (cz >= gz - 1) {
			cz = gz - 2;
		}
	} else {
		if (cx < 0) {
			cx = 0;
		}
		if (cy < 0) {
			cy = 0;
		}
		if (cz < 0) {
			cz = 0;
		}
		if (cx >= gx) {
			cx = gx - 1;
		}
		if (cy >= gy) {
			cy = gy - 1;
		}
		if (cz >= gz) {
			cz = gz - 1;
		}
	}
	
	 int idx = cx + cy * gx + cz * gxgy;
	// idx = ofClamp(idx, 0, numxyz-1);
	 
	return idx;	
}