#include "testApp.h"

float zspeed = 1.2;


//--------------------------------------------------------------
void testApp::setup(){
	
	ofSetFrameRate(60);		
	int mcres = 32;
	
	mc.setup(ofGetWidth(),ofGetHeight(),-ofGetWidth(), mcres, mcres, mcres);
		
//	ofEnableLighting();
	ofSetGlobalAmbientColor(ofColor(200,255,250));

	mc.isolevel = 0.02;
}

//--------------------------------------------------------------
void testApp::update(){
		

	if(ofGetMousePressed()){		
		float zpos = fmodf( ofGetFrameNum()*zspeed, mc.worlddim.z );
		zpos*=-1;
		ofPoint mp( ofGetMouseX(), ofGetMouseY(), zpos);
		ofPoint md( 100, 100, 100);
		mc.addCube(mp, md, 0.02);
		mc.polygoniseData();
	}
	
	
	
}

//--------------------------------------------------------------
void testApp::draw(){
//	cam.begin();
	ofBackground(0);
	ofSetColor(255);
	
	mc.draw();

	float zpos = fmodf( ofGetFrameNum()*-zspeed, mc.worlddim.z );
		
	ofPoint mp( ofGetMouseX(), ofGetMouseY(), zpos);
	
	if(ofGetMousePressed())
		ofSetColor(255,0,0);
	
	ofBox(mp.x, mp.y, mp.z, 20);
	
	ofSetColor(255,0,0);
	mc.drawnormals(0.021);
//	cam.end();
	
	ofSetColor(255);
	ofDrawBitmapString("fps: "+ofToString(ofGetFrameRate()) + "\n"+mc.getinfo()
					   +"\npos: "+ofToString(mp)
					   +"\n(c) clear"
					   +"\n(f) fade data"
					   +"\n(s) save stl",
					    10, 20);
}



//--------------------------------------------------------------
void testApp::keyPressed(int key){

	if(key=='s'){
		string fn = ofToString(mc.gx) +"x"+ofToString(mc.gy)+"x"+ofToString(mc.gz)+"-"+ofToString(ofGetFrameNum())+".stl";
		mc.saveStl(fn);
	}
	if(key=='c'){
		mc.clear();
		mc.polygoniseData();
	}
	if(key=='f'){
		mc.multData(0.9f);
		mc.polygoniseData();
	}
}

//--------------------------------------------------------------
void testApp::keyReleased(int key){

}

//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}