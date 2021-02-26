/*============================================================//
Copyright (c) 2020 America Makes
All rights reserved
Created under ALSAM project 3024

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//============================================================*/

/*============================================================//
ScanPath.cpp contains functions to generate scanpaths (hatches and 
contours) based on layer structure and configuration data.
Individual sets of scanpaths (corresponding to contours or
hatches from an individual part or region) make up Paths which are
appended to a vector of paths within a particular Trajectory.
A layer consists of one or more Trajectories which may each 
contain multiple Paths
//============================================================*/

#include "ScanPath.h"
#include "constants.h"

void findHatchBoundary(vector<vertex> &in, double hatchAngle, double *a_min, double *a_max)
{
	// Computes min/max x or y intersection of hatch lines drawn through all vertices
	// This identifies the number of hatch lines we need to create
	// Result is returned in a_min, a_max
	//
	// If hatchAngle 315 to 45 degrees or 135 to 225 degrees, we evaluate along the y axis
	// If hatchAngle is 45 to 135 or 225 to 315 degrees we evaluate along the x axis
	double hatchAngle_rads = hatchAngle * 3.14159265358979323846 / 180;  // hatch angle in radians
	int vLen = in.size();
	vector<double> a;
	a.reserve(vLen);  // pre-size the vector
	vertex v;

	// Determine whether to evaluate along x or y axis, based on hatchAngle
	if ((((int)hatchAngle + 315) % 180) > 90) {
		// hatches are spaced along y axis, so find y-intercept of hatchline through each vertex
		double slope = tan(hatchAngle_rads);
		//cout << "\nAlong y, hatchAngle " << hatchAngle << ", slope = " << slope << endl;
		for (int i = 0; i < vLen; i++)
		{
			v = in[i];
			a[i] = v.y - (v.x * slope);
			//cout << "   vx " << v.x << ", vy " << v.y << ", a[i] = " << a[i] << endl;  // testing purposes
		}
	}
	else {
		// hatches are spaced along x axis, so find x-intercept of hatchline through each vertex
		double slope = cos(hatchAngle_rads) / sin(hatchAngle_rads);  // can't use 1/tan, which fails at 90 degrees, whereas this returns 0
		//cout << "\nAlong x, hatchAngle " << hatchAngle << ", slope = " << slope << endl;
		for (int i = 0; i < vLen; i++)
		{
			v = in[i];
			a[i] = v.x - (v.y * slope);
			//cout << "   vx " << v.x << ", vy " << v.y << ", a[i] = " << a[i] << endl;  // testing purposes
		}
	}
	
	// find min/max values in a[] and return them in a_min, a_max
	*a_min = a[0];
	*a_max = a[0];
	for (int i = 0; i < vLen; i++)
	{	// using this rather than MinMax function due to some odd issues observed
		if (a[i] < *a_min) { *a_min = a[i]; }
		if (a[i] > *a_max) { *a_max = a[i]; }
	}
	//cout << "   a_min " << *a_min << ", a_max " << *a_max << endl;

	return;
}

double dist(vertex &v1, vertex &v2)
{
	double r = sqrt(pow((v1.x - v2.x), 2) + pow((v1.y - v2.y), 2));
	return r;
}

int findIntersection(vertex *out, double hatchAngle, double hatchAngle_rads, vector<vertex> BB, double intercept, edge e, double hatchFunctionValue)
{
	// e is the edge (2-point line segment) to be evaluated for intersection with a particular hatch line
	// the hatch line is defined via hatchAngle and intercept, which may be either x or y intercept based on hatchAngle
	// Return value:  0 indicates an intersection.  out will be populated with the intersect's x/y coordinates
	//				  Other values indicate no intersection.  out is undefined
	//				  
	// if hatchAngle is [-45 to 45] or [135 to 225] degrees, we evaluate along the y axis (intercept = y intercept)
	// otherwise, we evaluate along the x axis (intercept = x intercept)

	vertex edgeStart, edgeFinish, hatchStart, hatchFinish;
	double bL, bR, bB, bT;	//left, right, top, bottom coordinates of bounding box
	bL = BB[0].x;
	bR = BB[1].x;
	bB = BB[2].y;
	bT = BB[3].y;
	
	//vertices of the edge to be evaluated
	edgeStart = e.s;
	edgeFinish = e.f;
	
	//get vertices of the hatch segment to be evaluated based on hatch angle, intercept and the layer's bounding box
	// there are four possibilities:
	//	a) intercept may represent either x or y axis, depending on hatchAngle
	//	b) hatch direction may be positive or negative (along intercept angle), also depending on hatch Angle
	//
	// hatchAngle = 315 to 45 --> intercept is along the y axis and hatch runs in positive direction
	// hatchAngle = 135 to 225 = y axis, negative direction
	//
	// hatchAngle = 45 to 135 = intercept is along the x axis and hatch runs positive direction
	// hatchAngle = 225 to 315 = x axis, negative direction
	//
	if ( ( ( (int)hatchAngle + 315) % 180) > 90) {
		// hatchAngle is [315 to 45] or [135 to 225] degrees, i.e. closer to the x than y axis (and potentially parallel to the x axis)
		// therefore the hatches will be spaced along y axis, and intercept is assumed to be the y intercept of the current hatch line
		//
		// determine hatch direction
		if (hatchAngle > 90) {
			// angle is 135 to 225 degrees; hatch will run in negative y direction
			hatchStart.x = bR;
			hatchFinish.x = bL;
			hatchStart.y = intercept + hatchStart.x * hatchFunctionValue;		// y-coordinate on "left" side = y-intercept + tan(hatchAngle)*bL
			hatchFinish.y = intercept + hatchFinish.x * hatchFunctionValue;	// y-coordinate on "right" side = y-intercept + tan(hatchAngle)*bR
		}
		else {
			// hatch angle is 315 to 45 degrees; hatch will run in positive y direction
			hatchStart.x = bL;
			hatchFinish.x = bR;
			hatchStart.y = intercept + hatchStart.x * hatchFunctionValue;		// y-coordinate on "left" side = y-intercept + tan(hatchAngle)*bL
			hatchFinish.y = intercept + hatchFinish.x * hatchFunctionValue;	// y-coordinate on "right" side = y-intercept + tan(hatchAngle)*bR
		}
	}
	else {
		// hatchAngle is [45 to 135] or [225 to 315] degrees
		// therefore the hatches will be spaced along x axis, and intercept is assumed to be the x intercept of the current hatch line
		//
		// determine hatch direction
		if (hatchAngle > 180) {
			// angle is 225 to 315 degrees; hatches will run in negative x direction
			hatchStart.y = bB;
			hatchFinish.y = bT;
			hatchStart.x = intercept + hatchStart.y * hatchFunctionValue;
			hatchFinish.x = intercept + hatchFinish.y * hatchFunctionValue;
		}
		else {
			// angle is 45 to 135 degrees; hatches will run in positive x direction
			hatchStart.y = bT;
			hatchFinish.y = bB;
			hatchStart.x = intercept + hatchStart.y * hatchFunctionValue;
			hatchFinish.x = intercept + hatchFinish.y * hatchFunctionValue;
		}
	}

	//now that we have endpoints for the hatch line and the edge in question, check whether there is an intersection between them
	double a = hatchFinish.y - hatchStart.y;
	double b = hatchStart.x - hatchFinish.x;
	double c = a*(hatchStart.x) + b*(hatchStart.y);
	// Line CD represented as a2x + b2y = c2
	double a1 = edgeFinish.y - edgeStart.y;
	double b1 = edgeStart.x - edgeFinish.x;
	double c1 = a1*(edgeStart.x) + b1*(edgeStart.y);
	double det = a*b1 - a1*b;
	if (abs(det) / dist(edgeStart, edgeFinish) < minDeterminant) {
		// no intersection
		return 1;
	}
	else {
		// compute intersection point of the lines which correspond to segments AB/CD intersection
		double x = (b1*c - b*c1) / det;
		double y = (a*c1 - a1*c) / det;
		// determine whether intersection point falls on segment AB
		double minX, maxX, minY, maxY;
		minX = min(hatchStart.x, hatchFinish.x) - intersectRange;
		maxX = max(hatchStart.x, hatchFinish.x) + intersectRange;
		minY = min(hatchStart.y, hatchFinish.y) - intersectRange;
		maxY = max(hatchStart.y, hatchFinish.y) + intersectRange;
		if ((x < minX) | (x > maxX) | (y < minY) | (y > maxY)) {
			// intersection point does not fall on the hatch segment
			return 1;
		}
		else {
			// determine whether intersection point also falls on the edge segment
			minX = min(edgeStart.x, edgeFinish.x) - intersectRange;
			maxX = max(edgeStart.x, edgeFinish.x) + intersectRange;
			minY = min(edgeStart.y, edgeFinish.y) - intersectRange;
			maxY = max(edgeStart.y, edgeFinish.y) + intersectRange;
			if ((x < minX) | (x > maxX) | (y < minY) | (y > maxY)) {
				// intersection point does not fall on the edge segment
				return 1;
			}
			else {
				// hatch and edge intersect at a point that falls on both segments
				(*out).x = x;
				(*out).y = y;
				return 0;
			}
		}  // end if intersection falls on hatch and edge segments
	}  // end if det==0

	return 1;  // this return should never be reached; we should complete in one of the endpoints above
}

// sort vertexList in order of ascending y coordinate, and ascending x-coordinate if y's are equal
void yAsc(vector<vertex> &vertexList)
{
	int vLen = vertexList.size();
	for (int i = 0; i < vLen; i++)
	{
		for (int j = 0; j < vLen - 1; j++)
		{
			if (vertexList[j].y > vertexList[j + 1].y)
			{
				vertex tmp;
				tmp = vertexList[j + 1];
				vertexList[j + 1] = vertexList[j];
				vertexList[j] = tmp;
			}
			else if (vertexList[j].y == vertexList[j + 1].y)
			{
				if (vertexList[j].x > vertexList[j + 1].x)
				{
					vertex tmp;
					tmp = vertexList[j + 1];
					vertexList[j + 1] = vertexList[j];
					vertexList[j] = tmp;
				}
			}
		}
	}
	return;
}

void yDsc(vector<vertex> &vertexList)
{
	int vLen = (int)vertexList.size();
	for (int i = 0; i < vLen; i++)
	{
		for (int j = 0; j < vLen - 1; j++)
		{
			if (vertexList[j].y < vertexList[j + 1].y)
			{
				vertex tmp;
				tmp = vertexList[j + 1];
				vertexList[j + 1] = vertexList[j];
				vertexList[j] = tmp;
			}
			else if (vertexList[j].y == vertexList[j + 1].y)
			{
				if (vertexList[j].x < vertexList[j + 1].x)
				{
					vertex tmp;
					tmp = vertexList[j + 1];
					vertexList[j + 1] = vertexList[j];
					vertexList[j] = tmp;
				}
			}
		}
	}
	return;
}

void xAsc(vector<vertex> &vertexList)
{
	int vLen = vertexList.size();
	for (int i = 0; i < vLen; i++)
	{
		for (int j = 0; j < vLen - 1; j++)
		{
			if (vertexList[j].x > vertexList[j + 1].x)
			{
				vertex tmp;
				tmp = vertexList[j + 1];
				vertexList[j + 1] = vertexList[j];
				vertexList[j] = tmp;
			}
			else if (vertexList[j].x == vertexList[j + 1].x)
			{
				if (vertexList[j].y > vertexList[j + 1].y)
				{
					vertex tmp;
					tmp = vertexList[j + 1];
					vertexList[j + 1] = vertexList[j];
					vertexList[j] = tmp;
				}
			}
		}
	}
	return;
}

void xDsc(vector<vertex> &vertexList)
{
	int vLen = (int)vertexList.size();
	for (int i = 0; i < vLen; i++)
	{
		for (int j = 0; j < vLen - 1; j++)
		{
			if (vertexList[j].x < vertexList[j + 1].x)
			{
				vertex tmp;
				tmp = vertexList[j + 1];
				vertexList[j + 1] = vertexList[j];
				vertexList[j] = tmp;
			}
			else if (vertexList[j].x == vertexList[j + 1].x)
			{
				if (vertexList[j].y < vertexList[j + 1].y)
				{
					vertex tmp;
					tmp = vertexList[j + 1];
					vertexList[j + 1] = vertexList[j];
					vertexList[j] = tmp;
				}
			}
		}
	}
	return;
}

vector<vertex> eliminateDuplicateVertices(vector<vertex> &vertexList)
{	// assumes vertexList has been sorted in x and y.  doesn't matter which is the primary axis
	int vLen = (int)vertexList.size();
	vector<int> duplicateList;
	bool bContinue = true;
	int i = 0;

	// iterate through all vertices to identify those which are very close to the following vertex
	while (bContinue)
	{
		if ( (abs(vertexList[i].y - vertexList[i + 1].y) < overlapRange) & (abs(vertexList[i].x - vertexList[i + 1].x) < overlapRange))
		{	// vertices i and i+1 are identical within accuracy limit indicated by overlapRange.  save i to list for elimination
			duplicateList.push_back(i);
			i = i + 1;
		}
		else
			i = i + 1;
		if (i >= vLen - 1)
			bContinue = false;
	}

	// eliminate any duplicates and return as a new vector
	vector<vertex> vListOut = vertexList;
	for (int it = duplicateList.size() - 1; it >= 0; it--)
	{	// erase the duplicate vertex
		vListOut.erase(vListOut.begin() + duplicateList[it]);
	}
	return vListOut;
}

ray e2r(edge e)
{	
	ray r;
	r.x = (e.f).x - (e.s).x;
	r.y =(e.f).y-(e.s).y;
	double mod = rMod(r);
	r.x = r.x / mod;
	r.y = r.y / mod;
	return r;
}

ray rAdd(ray r1, ray r2)
{
	ray r;
	r.x = r1.x  + r2.x ;
	r.y = r1.y + r2.y;
	double mod = rMod(r);
	r.x = r.x / mod;
	r.y = r.y / mod;
	return r;
}

double rMod(ray r)
{
	return sqrt(pow(r.x, 2) + pow(r.y, 2));
}

double rAngle(ray r1, ray r2)
{
	double cTheta;
	cTheta = r1.x*r2.x + r1.y*r2.y;
	return acos(cTheta);
}

int getTurnDir(edge ev1, edge ev2)
{
	ray r1 = e2r(ev1);
	ray r2 = e2r(ev2);
	//take cross product to see if it is turning CW or CCW
	double p = r1.x*r2.y - r2.x*r1.y;
	if (p >= 0)
		return 1;
	else
		return -1;
}

void edgeOffset(layer &L, vector<int> regionIndex, vector<edge> &edgeListOut, vector<vector<edge>> &polyVectorsOut, double offset, bool returnPolyVectors)
{
	// this function offsets a set of edges "inward" for positive offset values
	// therefore outer contours will be indended toward center of part (making the outer contour smaller),
	// while inner contours (holes) will be indented in opposite fashion to make the hole larger
	//
	// uses the ClipperOffset routine, which may result in multiple output polygons if sections are clipped down to intersection
	//
	// regionIndex indicates the various regions (in layer L) to be included in the same offsetting operation.
	// L.s.rList[x].eList are the actual edges, where x is an element in regionIndex
	// **only one of edgeListOut or polyVectorsOut will be populated, based on returnPolyVectors:
	//		true --> polyVectorsOut is used, false --> edgeListOut instead
	// edgeListOut is a single vector which aggregates edges from all polygons that result from the operation
	// polyVectorsOut is a vector of edgeListOut-style edge vectors
	// offset is a value in millimeters.  positive shrinks the part (we flip the sign when passing it to ClipperOffset

	ClipperLib::Paths allContoursIn, allContoursOut;  // we expect mutiple contours (paths) in and out.  Counts may differ if some are eliminated by offsetting
	
	ClipperLib::Path contourIn;  // temporary aggregation variables
	ClipperLib::ClipperOffset coTemp;
	ClipperLib::IntPoint a;
	int i, j;	// iterators

	// ensure the output has no existing elements
	polyVectorsOut.clear();
	polyVectorsOut.shrink_to_fit();
	edgeListOut.clear();
	edgeListOut.shrink_to_fit();

	// 1. Iterate over regions and convert the edges into a Clipper Path (vector of endpoints).  Then add them to allContoursIn
	for (vector<int>::iterator it = regionIndex.begin(); it != regionIndex.end(); ++it) {
		// pre-size contourIn to hold the number of endpoints expected
		contourIn.clear();
		contourIn.shrink_to_fit();
		contourIn.reserve(L.s.rList[(*it)].eList.size());

		// iterate over the edges in this region and create a Clipper Path
		// Slic3r automatically reverses the order of holes (versus outer contours), so 
		// we can feed the vertices into ClipperLib directly, and they'll be interpreted accordingly
		for (i = 0; i < L.s.rList[(*it)].eList.size(); i++) {
			// push the starting point of this edge onto contourIn
			a.X = (int)round(L.s.rList[(*it)].eList[i].s.x / intersectRange);
			a.Y = (int)round(L.s.rList[(*it)].eList[i].s.y / intersectRange);
			contourIn.push_back(a);
		}
		// append contourIn to the vector of Paths
		allContoursIn.push_back(contourIn);
	}
	
	// 2. add path(s) to coTemp and set offsetting "styles" (but not the offset value itself)
	// appropriate settings for joinType (second parameter) are jtSquare or jtMiter.  jtRound is too fine in detail (adds many points)
	coTemp.AddPaths(allContoursIn, ClipperLib::jtMiter, ClipperLib::etClosedPolygon);

	// 3. offset allContoursIn --> allContoursOut.  This may lead to zero, one or multiple polygons
	coTemp.Execute(allContoursOut, -1*offset/intersectRange);  // must convert offset to same units as points, and flip sign because negative = shrinkage
	// allContoursOut is populated with the result, a vector of polygon-vertex vectors

	// 4. Determine if any polygons survived the offsetting operation
	if (allContoursOut.size() == 0) {
		// nothing left!
		return;
	}

	// 4a. we will begin by populating polyVectorsOut even if edgeListOut is desired instead, because we can create edgeListOut from polyVectorsOut
	// (but not the other way around).  Simplifies things
	polyVectorsOut.reserve(allContoursOut.size());  // assume all output polygons are valid (non-zero size)

	// 4b. iterate through allContoursOut. create a set of edges for each polygon and append it to polyVectorsOut
	vector<edge> tmpEdgeList;
	edge e;
	for (j=0; j < allContoursOut.size(); j++) {
		int polySize = allContoursOut[j].size();  // used multiple times
		if (polySize > 0) {
			tmpEdgeList.clear();	// clear any prior iterations
			tmpEdgeList.shrink_to_fit();
			tmpEdgeList.reserve(polySize);

			// iterate over this single polygon and create edges up to (but not including) final edge
			for (i = 0; i < polySize - 1; i++) {  // create edges connecting each point to its successor, except for final point
				e.s.x = allContoursOut[j][i].X * intersectRange;  // starting point of edge = this point
				e.s.y = allContoursOut[j][i].Y * intersectRange;

				e.f.x = allContoursOut[j][i + 1].X * intersectRange;  // ending point of edge = next point
				e.f.y = allContoursOut[j][i + 1].Y * intersectRange;

				tmpEdgeList.push_back(e);
			}
			// add the final edge, which connects the final point to starting point
			e.s.x = allContoursOut[j][polySize - 1].X * intersectRange;  // starting point of edge = final polygon point
			e.s.y = allContoursOut[j][polySize - 1].Y * intersectRange;
			e.f.x = allContoursOut[j][0].X * intersectRange;  // ending point of edge = first polygon point
			e.f.y = allContoursOut[j][0].Y * intersectRange;
			tmpEdgeList.push_back(e);

			// add tmpEdgeList polyVectorsOut as a new, separate, vector
			polyVectorsOut.push_back(tmpEdgeList);
		}  // if polySize > 0
	}  // end iteration over output polygons

	// 5. determine if the user wants the vector of vectors or simply an aggregation of all edges in one structure
	if (returnPolyVectors) {
		// no further work required.  just ensure the other return variable is cleared
		return;  // polyVectorsOut contains the result
	}
	else {
		// convert polyVectorsOut to a single set of edges
		// first, pre-size edgeList out based on the known number of edges
		int totalPts = 0;
		for (j=0; j < polyVectorsOut.size(); j++) { totalPts += polyVectorsOut[j].size(); }

		edgeListOut.reserve(totalPts);

		// iterate over polyVectorsOut (again) and append edges from each vector onto edgeListOut
		for (j=0; j < polyVectorsOut.size(); j++) {
			edgeListOut.insert(edgeListOut.end(), polyVectorsOut[j].begin(), polyVectorsOut[j].end());
		}
		return;  // edgeListOut contains the output edges
	}

	return; // should not reach this point
}

void displayPath(path P)
{
	/* To properly differentiate jumps and marks, need to reference the indicated segment style,
	find it in the list of styles and then determine if it either has NO traveler section, 
	or has a traveler with power=0
	*/
	vertex s, f;
	vector<segment> SP = P.vecSg;
	for (vector<segment>::iterator it = SP.begin(); it != SP.end(); ++it)
	{
		s = (*it).start;
		f = (*it).end;
/*		if ((*it).power == 0)  
			cout << "Jump: ";
		else
			cout << "Mark: ";  */
		cout << "Need to differentiate Jump and Mark: ";
		cout << s.x << " " << s.y << " ---> " << f.x << " " << f.y << endl;
	}
	cout << "Scan Files saved. " << endl;
	printf("====================================================================\n");
	return;
}

vector<vertex> getBB(layer &L)
{
	// If the layer is blank (contains no vertices), we return a default bounding box.
	// This might be the case when the layer consists only of single stripes and no STL parts
	//
	vector<vertex> BB;  // bounding box output; coordinates which define the four outermost x/y coordinates of the layer
	BB.reserve(4);
	vertex vTmp, vL, vR, vT, vB;
	vertex va, vb, vc, vd;
	vector<vertex>* vList = &(L.vList);

	if ((*vList).size() > 0) {
		
		vTmp.x = (*vList)[0].x;  // initialize the min/max vectors to the first point in vList
		vTmp.y = (*vList)[0].y;	 // if vList is empty, this will fail and go to the catch routine
		vL = vTmp;  // min x
		vR = vTmp;  // max x
		vT = vTmp;  // max y
		vB = vTmp;  // min y
		
		// iterate over every vertex in vList and compare to the left/right/top/bottom vars
		for (int k = 1; k < (*vList).size(); k++)
		{
			vTmp.x = (*vList)[k].x;
			vTmp.y = (*vList)[k].y;
			
			if (vTmp.x < vL.x)
				vL = vTmp;
			if (vTmp.x > vR.x)
				vR = vTmp;
			if (vTmp.y < vB.y)
				vB = vTmp;
			if (vTmp.y > vT.y)
				vT = vTmp;
		}
	} else {
		// Return a set of default bounding box dimensions
		vL.x = -10.0;
		vL.y = -10.0;
		vR.x = 10.0;
		vR.y = 10.0;
		vB.x = -10.0;
		vB.y = -10.0;
		vT.x = 10.0;
		vT.y = 10.0;
		//std::cout << "using default bounding box" << endl;
	}

	BB.push_back(vL);
	BB.push_back(vR);
	BB.push_back(vB);
	BB.push_back(vT);

	return BB;
}

int findInt(vector<vertex> &BB, vertex v1, vertex v2)
{
	int bInt = 0;

	if (v1.x < (BB[0].x - intersectRange) || v2.x < (BB[0].x - intersectRange))
		bInt = 1;
	if (v1.x > (BB[1].x + intersectRange) || v2.x > (BB[1].x + intersectRange))
		bInt = 1;
	if (v1.y < (BB[2].y - intersectRange) || v2.y < (BB[2].y - intersectRange))
		bInt = 1;
	if (v1.y > (BB[3].y + intersectRange) || v2.y > (BB[3].y + intersectRange))
		bInt = 1;
	return bInt;
}

vector<int> singleStripeCount(int layerNum, AMconfig &configData)
{
	// Evaluate the list of stripes to see if any are to be written on this layer and/or higher layers.
	// Return value is a vector of the trajectories of single stripes on this layer.
	// Knowing whether stripes exist on this layer allows us to predefine one or more stripe trajectories

	vector<int> stripeTrajThisLayer;
	int remainingStripeCount = 0;  // stripes remaining across this layer or higher
	int stripeCountThisLayer = 0;  // stripes on this layer only
	for (int stp = 0; stp < configData.stripeList.size(); stp++)
	{
		if ((configData.stripeList[stp].marked == false) & (configData.stripeList[stp].stripeLayerNum >= layerNum))
		{
			// stripe hasn't been marked in a prior layer. it remains for this or higher layer
			remainingStripeCount++;
			if (configData.stripeList[stp].stripeLayerNum == layerNum)
			{
				// stripe is to be marked on this layer
				stripeCountThisLayer++;
				// add its trajectory# to stripeTrajThisLayer (we'll sort and erase duplicates later)
				stripeTrajThisLayer.push_back(configData.stripeList[stp].trajectoryNum);
			}
		}
	}

	// if remainingStripeCount == 0, set configData.allStripesMarked = true and return
	if (remainingStripeCount == 0)
	{
		configData.allStripesMarked = true;
		return stripeTrajThisLayer;  // no stripes to mark on this layer (or any future layers)
	}

	// if all the remaining stripes are on this layer, set configData.allStripesMarked = true to avoid checking future layers
	if (remainingStripeCount <= stripeCountThisLayer) {
		configData.allStripesMarked = true;
	}

	// sort stripeTrajThisLayer in ascending order and remove duplicates
	sort(stripeTrajThisLayer.begin(), stripeTrajThisLayer.end());
	stripeTrajThisLayer.erase(unique(stripeTrajThisLayer.begin(), stripeTrajThisLayer.end()), stripeTrajThisLayer.end());

	return stripeTrajThisLayer;  // if we get to this point, stripeTrajThisLayer will be null if all remaining stripes are on higher layers; otherwise at least one element
}

path singleStripes(int layerNum, int trajectoryNum, AMconfig &configData)
{
	/* PROCESS:
		1. set stripesMarked = 0
		2. define an end coordinate vertex
		2. iterate over stripeList to find those which are: on this layer, in the desired trajectory, and not yet marked
				if stripesMarked > 0, add jump segment starting at prior segment end vertex and terminating at the start coordinate
				add mark segment from start/end of current stripe
				increment stripesMarked
	*/
	segment sg;				// holds a single segment for processing
	vector<segment> vSg;    // holds the output stripe & jump segments
	int stripesMarked = 0;  // running count of stripes marked on this layer
	vertex priorEndpoint;   // permits a jump from previous stripe, if any
	string jumpStyle;		// segment styles for all jumps between single stripes
	if (configData.outputIntegerIDs == true) {
		jumpStyle = to_string(configData.stripeJumpSegStyleIntID);
	}
	else {
		jumpStyle = configData.stripeJumpSegStyleID;
	}

	// Iterate across stripe list
	for (int stp = 0; stp < configData.stripeList.size(); stp++)
	{
		if ((configData.stripeList[stp].marked == false) & (configData.stripeList[stp].stripeLayerNum == layerNum) & (configData.stripeList[stp].trajectoryNum == trajectoryNum))
		{
			configData.stripeList[stp].marked = true;  // prevent this stripe from being marked again
			if (stripesMarked > 0)
			{
				// Add a jump from prior stripe to start of this one
				sg.start = priorEndpoint;
				sg.end.x = configData.stripeList[stp].startX;
				sg.end.y = configData.stripeList[stp].startY;
				sg.idSegStyl = jumpStyle;
				sg.isMark = 0;
				vSg.push_back(sg);
			}
			// Create a mark segment for the current stripe
			stripesMarked++;
			sg.start.x = configData.stripeList[stp].startX;
			sg.start.y = configData.stripeList[stp].startY;
			sg.end.x = configData.stripeList[stp].endX;
			sg.end.y = configData.stripeList[stp].endY;
			if (configData.outputIntegerIDs == true) {
				sg.idSegStyl = to_string(configData.stripeList[stp].segmentStyleIntID);
			}
			else {
				sg.idSegStyl = configData.stripeList[stp].segmentStyleID;
			}
			sg.isMark = 1;
			vSg.push_back(sg);
			priorEndpoint = sg.end;
		}
	}
	// Convert the segments into a path structure
	path P;
	P.vecSg = vSg;
	P.tag = configData.stripeRegionTag;
	P.type = "single_stripes";
	P.SkyWritingMode = configData.stripeSkywrtgMode;
	return P;
}

// CUSTOM ISLAND SCAN PATH CODE STARTS HERE

// Generates a randomized island order given the number of horizontal and vertical grid units (gridWidth and gridHeight respectively). 
// vertexOrder: A pointer to a vector of vertices. The value of each vertex is a certain island grid coordinate, and the order is randomized. We will iterate through this right after return.
// gridWidth: The number of grid units horiziontally.
// gridHeight: The number of grid units vertically.
void genIslandOrder(vector<vertex> *vertexOrder, int gridWidth, int gridHeight)
{
	// Create array with all possible data points
	for (int i = 0; i < gridWidth; i++) {
		for (int j = 0; j < gridHeight; j++) {
			vertex v;
			v.x = i;
			v.y = j;
			(*vertexOrder).push_back(v);
		}
	}

	// Run Fisher-Yates algorithm (Detailed @ https://www.geeksforgeeks.org/generate-a-random-permutation-of-1-to-n/) 
	// to randomize. 
	int indexToSwitch;
	vertex swap;
	for (int i = (*vertexOrder).size() - 1; i >= 0; i--) {

		// We want it to be [0, i] so we use rand() % (i + 1)
		indexToSwitch = rand() % (i + 1);

		// Swap it with the index we are switching with 
		swap = (*vertexOrder)[indexToSwitch];
		(*vertexOrder)[indexToSwitch] = (*vertexOrder)[i];
		(*vertexOrder)[i] = swap;
	}

	// No explicit return value; vertexOrder is effectively our return but we pass by pointer
}

// Stores bounding box as a TopLeft, TopRight, BottomLeft, BottomRight representation
// instead of their default array [Leftmost, Rightmost, Bottommost, Topmost]
struct BoundingBox {
	vertex tl, tr, bl, br;
};

// Got really annoying to work with all the points. Simpler representation of a segment.
struct SimpleSegment {
	vertex v1, v2;
};

// Deep copy vertex (used to avoid aliasing issues)
vertex copyVertex(vertex src) {
	vertex v;
	v.x = src.x;
	v.y = src.y;
	return v;
}

// Deep copy segment (used to avoid aliasing issues)
SimpleSegment copySegment(SimpleSegment src) {
	SimpleSegment segment;
	segment.v1.x = src.v1.x;
	segment.v1.y = src.v1.y;
	segment.v2.x = src.v2.x;
	segment.v2.y = src.v2.y;
	return segment;
}

/* Basically a toString() for all the data types we have. Useful for debugging. */
void printVertex(vertex v) { 
	if (v.x == NULL || v.y == NULL) {
		cout << "printVertex: Invalid vertex." << endl;
		return;
	}
	cout << "(" << v.x << "," << v.y << ")"; 
}
void printSegment(SimpleSegment segment) { cout << "((" << segment.v1.x << "," << segment.v1.y << "),(" << segment.v2.x << "," << segment.v2.y << "))"; }
void printBB(BoundingBox bb) { 
	cout << "(tl,tr,bl,br): (";
	printVertex(bb.tl); cout << ",";
	printVertex(bb.tr); cout << ",";
	printVertex(bb.bl); cout << ",";
	printVertex(bb.br); cout << ")";
}

inline double Det(double a, double b, double c, double d)
{
	return a*d - b*c;
}

/// Calculate intersection of two lines.
///\return true if found, false if not found or error
bool lineIntersectsLine(SimpleSegment *s1, SimpleSegment *s2, vertex *intersect)
{
	// Found this implementation online; easier to just translate my system of segments here than replace every one and risk screwing it up
	double x1 = s1->v1.x, y1 = s1->v1.y;
	double x2 = s1->v2.x, y2 = s1->v2.y;
	double x3 = s2->v1.x, y3 = s2->v1.y;
	double x4 = s2->v2.x, y4 = s2->v2.y;
	
	
	double detL1 = Det(x1, y1, x2, y2);
	double detL2 = Det(x3, y3, x4, y4);
	double x1mx2 = x1 - x2;
	double x3mx4 = x3 - x4;
	double y1my2 = y1 - y2;
	double y3my4 = y3 - y4;

	double denom = Det(x1mx2, y1my2, x3mx4, y3my4);
	if (denom == 0.0) // Lines don't seem to cross
	{
		intersect->x = NAN;
		intersect->y = NAN;
		return false;
	}

	double xnom = Det(detL1, x1mx2, detL2, x3mx4);
	double ynom = Det(detL1, y1my2, detL2, y3my4);
	intersect->x = xnom / denom;
	intersect->y = ynom / denom;
	if (!isfinite(intersect->x) || !isfinite(intersect->y)) // Probably a numerical issue
		return false;

	return true; //All OK
}

// Given their edge-limited version of boundingBox, converts it to my corner-based representation.
void getCustomBoundingBox(BoundingBox *ourBoundingBox, vector<vertex> *theirBoundingBox)
{
	(*ourBoundingBox).tl.x = (*theirBoundingBox)[0].x;
	(*ourBoundingBox).tl.y = (*theirBoundingBox)[3].y;
	(*ourBoundingBox).br.x = (*theirBoundingBox)[1].x;
	(*ourBoundingBox).br.y = (*theirBoundingBox)[2].y;
	(*ourBoundingBox).bl.x = (*ourBoundingBox).tl.x;
	(*ourBoundingBox).bl.y = (*ourBoundingBox).br.y;
	(*ourBoundingBox).tr.x = (*ourBoundingBox).br.x;
	(*ourBoundingBox).tr.y = (*ourBoundingBox).tl.y;
}

// Fills islandBoundingBox with the correct values given the overall bounding box, the specified island, and the overall island size.
// This becomes a little more complicated when we add the island offset.
void fillIslandCoordsAbsolute(BoundingBox *islandBoundingBox, BoundingBox *overallBoundingBox, vertex *island, float islandSize)
{
	(*islandBoundingBox).bl.x = (*overallBoundingBox).bl.x + (*island).x * islandSize;
	(*islandBoundingBox).bl.y = (*overallBoundingBox).bl.y + (*island).y * islandSize;
	(*islandBoundingBox).tl.x = (*islandBoundingBox).bl.x;
	(*islandBoundingBox).tl.y = (*islandBoundingBox).bl.y + islandSize;
	(*islandBoundingBox).br.x = (*islandBoundingBox).bl.x + islandSize;
	(*islandBoundingBox).br.y = (*islandBoundingBox).bl.y;
	(*islandBoundingBox).tr.x = (*islandBoundingBox).bl.x + islandSize;
	(*islandBoundingBox).tr.y = (*islandBoundingBox).br.y + islandSize;
}

// Creates a bounding box that has padAmount added onto each edge of the input bounding box
BoundingBox padBoundingBox(BoundingBox *inputBoundingBox, float padAmount) {
	BoundingBox output;
	output.tl.x = (*inputBoundingBox).tl.x - padAmount;
	output.tl.y = (*inputBoundingBox).tl.y + padAmount;
	output.tr.x = (*inputBoundingBox).tr.x + padAmount;
	output.tr.y = (*inputBoundingBox).tr.y + padAmount;
	output.bl.x = (*inputBoundingBox).bl.x - padAmount;
	output.bl.y = (*inputBoundingBox).bl.y - padAmount;
	output.br.x = (*inputBoundingBox).br.x + padAmount;
	output.br.y = (*inputBoundingBox).br.y - padAmount;
	return output;
}

bool verticesEqual(vertex v1, vertex v2) {
	return v1.x == v2.x && v1.y == v2.y;
}

bool pointInBB(vertex v1, BoundingBox bb) {
	return v1.x >= bb.tl.x && v1.x <= bb.tr.x && v1.y >= bb.bl.y && v1.y <= bb.tl.y;
}

SimpleSegment makeSegment(double x1, double y1, double x2, double y2) {
	SimpleSegment segment;
	segment.v1.x = x1;
	segment.v1.y = y1;
	segment.v2.x = x2;
	segment.v2.y = y2;
	return segment;
}

void printVector(vector<vertex> *v) {
	int i;
	for (i = 0; i < v->size(); i++) {
		printVertex((*v)[i]);
	}
}

// Used to slightly offset input point to make sure 
void fuzzPoint(vertex *v)
{
	float fuzzAmount = .01 * (((float)rand() / RAND_MAX) * 2 - 1);
	v->x += fuzzAmount;
	v->y += fuzzAmount;
}

// Given the overall vector and a given island, returns the subsegment of that vector that is inside the island. 
// THIS IS CALLED WHEN THERE IS NO INTERSECTION. If this is the case, the function will return false.
// Returns true if there was an intersecting vector, false if there wasn't.
bool getIntersectingVector(SimpleSegment *subsegment, SimpleSegment *segment, BoundingBox *islandBB)
{
	// Default is returning two vertices with null values.
	subsegment->v1.x = NULL;
	subsegment->v1.y = NULL;
	subsegment->v2.x = NULL;
	subsegment->v2.y = NULL;

	// Construct segments for each side of our island (for comparison purposes)
	SimpleSegment top, left, right, bottom;
	top.v1 = islandBB->tl; top.v2 = islandBB->tr;
	left.v1 = islandBB->tl; left.v2 = islandBB->bl;
	right.v1 = islandBB->tr; right.v2 = islandBB->br;
	bottom.v1 = islandBB->bl; bottom.v2 = islandBB->br;

	// Edge Case: Colinear to one of the edges. The function implementations I am using are a little undefined here so I handle them manually to be sure.
	// This will likely never happen due to floating point precision, but gotta be ready incase it does! 

	// IF segment is straight horizontal AND segment is straight horizontal to one of the vertical portions
	// Note that we don't consider anything right colinear, only left, to avoid duplicates entries across islands.
	if (segment->v1.y == segment->v2.y && segment->v1.y == islandBB->bl.y) {

		// If there is no overlap, get out.
		if ((segment->v1.x <= islandBB->tl.x && segment->v2.x <= islandBB->tl.x) ||
			segment->v1.x >= islandBB->tl.x && segment->v2.x >= islandBB->tl.x) {
			return false;
		}
		
		vertex leftPoint = segment->v1.x <= segment->v2.x ? segment->v1 : segment->v2;
		vertex rightPoint = segment->v1.x <= segment->v2.x ? segment->v2 : segment->v1;

		// Figure out left endpoint 
		double left;
		if (leftPoint.x >= islandBB->tl.x) {
			left = leftPoint.x;
		} else {
			left = islandBB->tl.x;
		}

		// Determine right endpoint
		double right;
		if (rightPoint.x <= islandBB->tr.x) {
			right = rightPoint.x;
		} else {
			right = islandBB->tr.x;
		}

		// Subsegment is the y we already established, then left to right boundaries
		*subsegment = makeSegment(left, segment->v1.y, right, segment->v1.y);
		return true;
	}

	// Specifically want to exclude top ones as previously mentioned (to avoid duplicates)
	if (segment->v1.y == segment->v2.y && segment->v1.y == islandBB->tl.y) {
		return false;
	}

	// Vertical Cases
	// Note that we don't consider anything top colinear, only bottom, to avoid duplicate entries across islands.
	if (segment->v1.x == segment->v2.x && segment->v1.x == islandBB->tl.x) {

		// If there is no overlap, get out.
		if ((segment->v1.y <= islandBB->bl.y && segment->v2.y <= islandBB->bl.y) ||
			segment->v1.y >= islandBB->tl.y && segment->v2.y >= islandBB->tl.y) {
			return false;
		}

		vertex bottomPoint = segment->v1.y <= segment->v2.y ? segment->v1 : segment->v2;
		vertex topPoint = segment->v1.y <= segment->v2.y ? segment->v2 : segment->v1;

		// Figure out left endpoint 
		double bottom;
		if (bottomPoint.y >= islandBB->bl.y) {
			bottom = bottomPoint.y;
		} else {
			bottom = islandBB->bl.y;
		}

		// Determine right endpoint
		double top;
		if (topPoint.y <= islandBB->tl.y) {
			top = topPoint.y;
		} else {
			top = islandBB->tl.y;
		}

		// Subsegment is the y we already established, then left to right boundaries
		*subsegment = makeSegment(segment->v1.x, bottom, segment->v1.x, top);
		return true;
	}

	// Specifically want to exclude right ones as previously mentioned (to avoid duplicates)
	if (segment->v1.x == segment->v2.x && segment->v1.x == islandBB->tr.x) {
		return false;
	}

	// Check each side. If there is an intersection we know the intersection point will be filled, so we can then check that it's actually within the segment.
	// Again, note that lineIntersectsLine just returns if they *ever* intersect, and if they do it returns a point. If we know they intersect we have to compound it with checking that it's within our bounds.
	vertex topIntersection, rightIntersection, bottomIntersection, leftIntersection;
	bool intersectsTop = lineIntersectsLine(&top, segment, &topIntersection) && // They aren't parallel
		(topIntersection.x >= islandBB->tl.x && topIntersection.x <= islandBB->tr.x) && // Their intersection point is on the actual island boundary itself 
		((segment->v1.y <= islandBB->tl.y && segment->v2.y >= islandBB->tl.y) || // This and last condition ensure segment is actually split across, not just intersecting if we extend it
		(segment->v1.y >= islandBB->tl.y && segment->v2.y <= islandBB->tl.y));
	bool intersectsRight = lineIntersectsLine(&right, segment, &rightIntersection) &&
		(rightIntersection.y >= islandBB->bl.y && rightIntersection.y <= islandBB->tl.y) &&
		((segment->v1.x <= islandBB->tr.x && segment->v2.x >= islandBB->tr.x) ||
		(segment->v1.x >= islandBB->tr.x && segment->v2.x <= islandBB->tr.x));
	bool intersectsBottom = lineIntersectsLine(&bottom, segment, &bottomIntersection) && 
		(bottomIntersection.x >= islandBB->tl.x && bottomIntersection.x <= islandBB->tr.x) &&
		((segment->v1.y <= islandBB->bl.y && segment->v2.y >= islandBB->bl.y) || 
		(segment->v1.y >= islandBB->bl.y && segment->v2.y <= islandBB->bl.y));
	bool intersectsLeft = lineIntersectsLine(&left, segment, &leftIntersection) &&
		(leftIntersection.y >= islandBB->bl.y && leftIntersection.y <= islandBB->tl.y) && 
		((segment->v1.x <= islandBB->tl.x && segment->v2.x >= islandBB->tl.x) ||
		(segment->v1.x >= islandBB->tl.x && segment->v2.x <= islandBB->tl.x));

	int numIntersections = 0;
	if (intersectsTop) numIntersections++; 
	if (intersectsRight) numIntersections++; 
	if (intersectsBottom) numIntersections++; 
	if (intersectsLeft) numIntersections++; 

	// If no intersections
	if (numIntersections == 0) {

		// It's either both are in or neither is in at this point. If the first is out, whole segment is out.
		if (segment->v1.x < islandBB->tl.x || segment->v1.x > islandBB->tr.x ||
			segment->v1.y < islandBB->bl.y || segment->v1.y > islandBB->tl.y) {
			return false;
		}

		// Otherwise, whole segment is in.
		*subsegment = copySegment(*segment); // Normally I'd just switch the pointer but we want to copy it all over to avoid aliasing 
	}

	// If one intersection, the line is the one intersection point we know plus whichever endpoint is within the rectangle.
	else if (numIntersections == 1) {

		// Figure out which intersection is inside the rectangle
		if (segment->v1.x >= islandBB->tl.x && segment->v1.x <= islandBB->tr.x && 
			segment->v1.y >= islandBB->bl.y && segment->v1.y <= islandBB->tl.y) {
			subsegment->v1 = copyVertex(segment->v1);
		}
		else {
			subsegment->v1 = copyVertex(segment->v2);
		}

		// Other half of the vector is from our original intersection calculations
		if (intersectsTop) {
			subsegment->v2 = copyVertex(topIntersection);
		}
		else if (intersectsRight) {
			subsegment->v2 = copyVertex(rightIntersection);
		}
		else if (intersectsBottom) {
			subsegment->v2 = copyVertex(bottomIntersection);
		}
		else if (intersectsLeft) {
			subsegment->v2 = copyVertex(leftIntersection);
		}
	}

	// If two intersections, it must intersect two sides of our rectangle. We add any of the two vertices that aren't duplicates.
	// Also covers 3 and 4 intersection cases because that just means a vertex double counts. We just make sure we pick two vertices of our set that aren't the same.
	else if (numIntersections >= 2 && numIntersections <= 4) {

		// Edge Case: One point on island corner and other point outside of island in a direction that has no overlap w/ island
		if (numIntersections == 2) {

			// Top-Right: Error out in case that intersection is exactly on and other point is outside island
			if (intersectsTop && intersectsRight &&
				((verticesEqual(segment->v1, islandBB->tr) && !pointInBB(segment->v2, *islandBB)) ||
				(verticesEqual(segment->v2, islandBB->tr) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check top-left
			if (intersectsTop && intersectsLeft &&
				((verticesEqual(segment->v1, islandBB->tl) && !pointInBB(segment->v2, *islandBB)) ||
					(verticesEqual(segment->v2, islandBB->tl) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check bottom-right
			if (intersectsBottom && intersectsRight &&
				((verticesEqual(segment->v1, islandBB->br) && !pointInBB(segment->v2, *islandBB)) ||
					(verticesEqual(segment->v2, islandBB->br) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}

			// Check bottom-left 
			if (intersectsBottom && intersectsLeft &&
				((verticesEqual(segment->v1, islandBB->bl) && !pointInBB(segment->v2, *islandBB)) ||
					(verticesEqual(segment->v2, islandBB->bl) && !pointInBB(segment->v1, *islandBB)))) {
				return false;
			}
		}

		bool firstFilled = false, secondFilled = false;

		// Handle top intersection 
		if (intersectsTop) {
			subsegment->v1 = copyVertex(topIntersection);
			firstFilled = true;
		}

		// Handle bottom intersection 
		if (intersectsBottom) {

			// See if first AND check that it's not a duplicate. 
			// Note that short-circuiting results in this never accessing something that hasn't been assigned its value in our system yet.
			if (firstFilled && !verticesEqual(subsegment->v1, bottomIntersection)) {
				subsegment->v2 = copyVertex(bottomIntersection);
				secondFilled = true;
			}

			else if (!firstFilled) {
				subsegment->v1 = copyVertex(bottomIntersection);
				firstFilled = true;
			}
		}

		// Handle right intersection
		if (!secondFilled && intersectsRight) {

			// See if second AND check that it's not a duplicate
			if (firstFilled && !verticesEqual(subsegment->v1, rightIntersection)) {
				subsegment->v2 = copyVertex(rightIntersection);
				secondFilled = true;
			}

			// Otherwise, see if first
			else if (!firstFilled)
			{
				subsegment->v1 = copyVertex(rightIntersection);
				firstFilled = true;
			}
		}

		// Handle left intersection
		if (!secondFilled && intersectsLeft) {

			// See if second AND check that it's not a duplicate
			if (firstFilled && !verticesEqual(subsegment->v1, leftIntersection)) {
				subsegment->v2 = copyVertex(leftIntersection);
				secondFilled = true;
			}

			// Otherwise, see if first
			else if (!firstFilled)
			{
				subsegment->v1 = copyVertex(leftIntersection);
				firstFilled = true;
			}
		}

		// Error case that provides more debug information
		if (!secondFilled) {
			cout << "getIntersectingVector: Got through 2/3/4 intersection code without two intersections!" << endl;
			cout << "numIntersections: " << numIntersections << endl;
			cout << "topIntersection: "; printVertex(topIntersection); cout << endl;
			cout << "bottomIntersection: "; printVertex(bottomIntersection); cout << endl;
			cout << "leftIntersection: "; printVertex(leftIntersection); cout << endl;
			cout << "rightIntersection: "; printVertex(rightIntersection); cout << endl;
			cout << "segment: "; printSegment(*segment); cout << endl;
			cout << "bb: "; printBB(*islandBB); cout << endl;
			return false;
		}
	}

	// Error case that provides more debug information
	else {
		cout << "getIntersectingVector: numIntersections was " << numIntersections << " instead of 0/1/2/3!" << endl;
		return false;
	}	

	// Assume we found one; We return false elsewhere if we hit any of the error conditions
	return true;
}

// Splits the given list of intersections based on bounding box size and island size
// Returns 0 if successful, 1 if unsuccessful
void splitByIslands(vector<vertex> *splitVectors, vector<vertex> *originalVectors, vector<vertex> *boundingBox, double islandSize, double islandOffset, double islandOverlap, double hatchAngle, int layerNum)
{
	// isList: List of vertices representing where the initial hatch lines intersect the part.
	// boundingBox: Contains the points representing the min/max x value and min/max y value. 
	// islandSize: The size we want our islands to be. If this is bigger than the part itself, then we just have one big island.

	// Convert their vectors to SimpleSegment, a simple data structure I use so I'm not always dealing with pairs of vertices
	// cout << "originalVectors (size " << originalVectors->size() << "): " << endl;
	vector<SimpleSegment> originalSegments;
	for (int i = 0; i < originalVectors->size(); i += 2) {
		SimpleSegment segment;
		segment.v1 = copyVertex((*originalVectors)[i]);
		segment.v2 = copyVertex((*originalVectors)[i + 1]);
		originalSegments.push_back(segment);
	}

	// Convert their bounding box representation to a corner-based representation, which works better for us
	BoundingBox overallBoundingBox;
	getCustomBoundingBox(&overallBoundingBox, boundingBox);

	// The easiest way to get our offset is to shift the bottom-left corner of our bounding box by the offset, 
	// then keep going up-right, stopping one before the first that is fully inside our original bbox.
	BoundingBox overallBoundingBoxCopy;
	overallBoundingBoxCopy.tl = copyVertex(overallBoundingBox.tl);
	overallBoundingBoxCopy.tr = copyVertex(overallBoundingBox.tr);
	overallBoundingBoxCopy.bl = copyVertex(overallBoundingBox.bl);
	overallBoundingBoxCopy.br = copyVertex(overallBoundingBox.br);

	// Apply offset, then step back up and to the right until we get to the smallest new bl without being smaller than original
	overallBoundingBox.bl.x -= islandOffset;
	overallBoundingBox.bl.y -= islandOffset;
	bool done = false;
	vertex cache = copyVertex(overallBoundingBox.bl);
	while (!done) {
		overallBoundingBox.bl.x += islandSize;
		overallBoundingBox.bl.y += islandSize;
		if (overallBoundingBox.bl.x >= overallBoundingBoxCopy.bl.x &&
			overallBoundingBox.bl.y >= overallBoundingBoxCopy.bl.y) {
			overallBoundingBox.bl = copyVertex(cache);
			done = true;
		}
		cache = copyVertex(overallBoundingBox.bl);
	}
	overallBoundingBox.tl.x = overallBoundingBox.bl.x;
	overallBoundingBox.br.y = overallBoundingBox.bl.y;

	// Various information commonly used throughout function 
	const int gridWidth = ceil((overallBoundingBox.tr.x - overallBoundingBox.tl.x) / islandSize); // Round up ( (Rightmost.x - Leftmost.x) / islandSize )
	const int gridHeight = ceil((overallBoundingBox.tr.y - overallBoundingBox.br.y) / islandSize); // Round up ( (Topmost.y - Bottommost.y) / islandSize )

	// Get a randomized order that we will go through the islands in 
	vector<vertex> islandOrder;
	genIslandOrder(&islandOrder, gridWidth, gridHeight);

	// Iterate through the island order we just generated 
	for (int i = 0; i < islandOrder.size(); i++) {

		// Ignore every "other" one (makes it so adjacent islands have perpendicular islands)
		// Because we call hatch twice, we have to do this based on island coordinates rather than iteration order (which will be random, and different through both hatch() calls)
		if (layerNum % 2 == 0) {
			if (hatchAngle == 0) {
				if (fmod(islandOrder[i].x, 2) == 0 && fmod(islandOrder[i].y, 2) == 0 ||
					fmod(islandOrder[i].x, 2) == 1 && fmod(islandOrder[i].y, 2) == 1) {
					continue;
				}
			} else if (hatchAngle == 90) {
				if (fmod(islandOrder[i].x, 2) == 0 && fmod(islandOrder[i].y, 2) == 1 ||
					fmod(islandOrder[i].x, 2) == 1 && fmod(islandOrder[i].y, 2) == 0) {
					continue;
				}
			} else {
				cout << "ERROR: Hatch angle was neither 0 or 90!" << endl;
			}
		} else {
			if (hatchAngle == 0) {
				if (fmod(islandOrder[i].x, 2) == 0 && fmod(islandOrder[i].y, 2) == 1 ||
					fmod(islandOrder[i].x, 2) == 1 && fmod(islandOrder[i].y, 2) == 0) {
					continue;
				}
			} else if (hatchAngle == 90) {
				if (fmod(islandOrder[i].x, 2) == 0 && fmod(islandOrder[i].y, 2) == 0 ||
					fmod(islandOrder[i].x, 2) == 1 && fmod(islandOrder[i].y, 2) == 1) {
					continue;
				}
			} else {
				cout << "ERROR: Hatch angle was neither 0 or 90!" << endl;
			}
		}

		// Get boundaries of this specific island 
		BoundingBox islandBB;
		fillIslandCoordsAbsolute(&islandBB, &overallBoundingBox, &(islandOrder[i]), islandSize);

		// Padded version is what we pass in as islandBB
		// BoundingBox paddedBBox = padBoundingBox(&(islandBB), 0);
		BoundingBox paddedBBox = padBoundingBox(&(islandBB), islandOverlap / 2);

		// Iterate through all segments. Add the subsegment that intersects to our end array.
		for (int j = 0; j < originalSegments.size(); j++) {
			SimpleSegment intersectingSegment;
			if (getIntersectingVector(&intersectingSegment, &(originalSegments[j]), &(paddedBBox)))
			{
				(*splitVectors).push_back(intersectingSegment.v1);
				(*splitVectors).push_back(intersectingSegment.v2);
			}
		}
	}
}

// CUSTOM ISLAND SCAN PATH CODE ENDS HERE

path hatch(layer &L, vector<int> regionIndex, regionProfile &rProfile, bool useIslandScanning, double islandSize, double islandOffset, double islandOverlap, double offset, double hatchAngle, double a_min, double a_max, bool outputIntegerIDs, vector<vertex> boundingBox, int layerNum)
{
	/*	L: pointer to the layer structure
	regionIndex: list of regions to be hatched.  Must all have the same tag identified in rProfile
		Should indicate all related regions *and* holes, because all outer edges will be assumed to contain surface area and not holes
	rProfile: pointer to the region profile for these regions
	offset: hatch offset (from rProfile) with all contour offsets added
	hatchAngle: hatch angle in degrees.  0 = horizontal to the left
	a_min, a_max: min/max intercepts of hatch lines with vertices in L.vList, based on whether hatchAngle is closer to x or y axis
		if hatchAngle is [315 to 45] or [135 to 225], these should be min/max y coordinates
		however if hatchAngle is [45 to 135] or [225 to 315], these should be min/max x coordinates
	outputIntegerIDs: whether to use auto-generated integer ID's for segments styles, or the original string ID's
	boundingBox: vector of min/max x and y coordinates of this layer
	*/

	path P;  // path of hatch segments to be created
	vector<edge> edgeList;
	vector<segment> vSg;

	double hSpace = rProfile.resHatch;  // hSpace is the spacing between hatches, which will be modified by hatchAngle
	double a_curr;  // starting intercept coordinate of the current hatch line

	int dirHatch = 0;  // controls hatch direction; 0 is positive hatch direction, 1 = negative hatch direction
	// current code flips dirHatch on every hatch line so that the hatches zig-zag.
	// OPTIONAL:  allow flipping of dirHatch to be controlled by user.  If disabled, this would cause hatches to all be drawn in one direction
	// OPTIONAL:  in addition, all user to set the value of dirHatch.  Not so useful, since dirHatch==1 would change behavior with hatch angle

	// Determine primary hatching orientation (toward x or y axis), and adjust hSpace to account for hatchAngle
	double hatchAngle_rads = hatchAngle * 3.14159265358979323846 / 180;  // hatch angle in radians
	double hatchFunctionValue;  // sine or tangent value of hatchAngle_rads to be computed once and passed to findIntersection
	string primaryHatchDir;     // either y (progress along y axis) or x (progress along x axis), depending on hatchAngle
	// Determine whether to generate hatches progressing along the vertical or horizontal axis based on hatch angle
	if ((((int)hatchAngle + 315) % 180) > 90) {
		// hatchAngle is [315 to 45] or [135 to 225] degrees, i.e. closer to the x than y axis (and potentially parallel to the x axis)
		// we will increment hatches along the y axis, which will never be parallel to the hatch direction
		primaryHatchDir = "y";
		hSpace = hSpace / cos(hatchAngle_rads);  // this is the hatch-to-hatch spacing upon intersection with y axis
		hatchFunctionValue = tan(hatchAngle_rads);
	}
	else {
		// hatchAngle is [45 to 135] or [225 to 315] degrees, i.e. closer to the y than x axis (and potentially parallel to the y axis)
		// we will increment hatches along the x axis, which will never be parallel to the hatch direction
		primaryHatchDir = "x";
		hSpace = hSpace / sin(hatchAngle_rads);  // this is the hatch-to-hatch spacing upon intersection with x axis
		hatchFunctionValue = cos(hatchAngle_rads) / sin(hatchAngle_rads);  // can't use 1/tan, which fails at 90 degrees;
	}

	int numRegions = regionIndex.size();
	int naiveTmpSize = numRegions * 3;  // estimated number of intersections per hatch line, for sizing temporary vectors

	// Identify the segment styles, in terms of either their string or auto-generated integer ID's
	string hatchSegStyle, jumpSegStyle;
	if (outputIntegerIDs == true) {
		hatchSegStyle = to_string(rProfile.hatchStyleIntID);
		jumpSegStyle  = to_string(rProfile.jumpStyleIntID);
	}
	else {
		hatchSegStyle = rProfile.hatchStyleID;
		jumpSegStyle  = rProfile.jumpStyleID;
	}

	//*** Offset the edges across all regions under this tag (per hatch offset) and combine into a list of offsetted edges across all regions using this tag
	vector<vector<edge>> unusedVar;  // alternate output for edgeOffset, which must exist even if it won't be populated
	edgeOffset(L, regionIndex, edgeList, unusedVar, offset, false);  // false=return all offset edges in one (jumbled) vector

	// determine whether we have anything to output (any polygons which survived offsetting)
	if (edgeList.size() == 0) {
		// no result!
		P.type = "";  // indicates null result
		P.vecSg.clear();
		return P;
	}

	// *** Identify intersection points between edges and parallel-spaced hatch lines by stepping from a_min to a_max in increments of hSpace.
	double a_start, a_end;
	if (hSpace > 0) {
		a_start = a_min;
		a_end = a_max;
	}
	else {
		a_start = a_max;
		a_end = a_min;
	}

	bool finished = false;  // use this to determine when we've reached a_end, since we could be incrementing or decrementing via hSpace
	// We start at a_start + hSpace because a hatch exactly at a_start would have no length (it would be a one-point hatch)
	vector<vertex> isList;  // list of intersections between edges and hatch lines
	a_curr = a_start + hSpace;  // a_curr is the x or y intersection of the current hatch line
	while (!finished)
	{
		vector<vertex> tmp_isList, tmp_isListNoDuplicates;  // temporary list of intersections between one hatch line and all edges for this region tag
		tmp_isList.reserve(naiveTmpSize);  // number of intersections is not yet known; may be 0 or very large depending on geometry
		for (vector<edge>::iterator it = edgeList.begin(); it != edgeList.end(); ++it)
		{
			//get intersection with edge
			vertex is;
			int isB;
			isB = findIntersection(&is, hatchAngle, hatchAngle_rads, boundingBox, a_curr, *it, hatchFunctionValue);

			// if an intersection was found (return value == 0), add to tmp intersection vector
			if (isB == 0)
				tmp_isList.push_back(is);
		}

		// Now that we've evaluated a particular hatch line against all edges, process whatever intersections we have
		//
		// If there are intersections, sort them, then eliminate any duplicates
		if (tmp_isList.size() > 0) {
			// Sort tmp_isList in x or y-direction, according to hatch direction
			//
			// Note - this method will always begin with intersection with an outer edge and proceed across
			// an area to be marked.  No matter how many intersections are crossed, the tmp_isList
			// will alternate between mark and jump segments.
			// For instance, given six vertices in tmp_isList (0-5), the resulting edges will form
			// segments 0-1 (mark), 1-2 (jump), 2-3 (mark), 3-4 (jump), 4-5 (mark), followed by a jump
			// to the next vertex, if any
			if (dirHatch == 0)
			{
				dirHatch = 1;
				if (primaryHatchDir == "x") {
					yAsc(tmp_isList);  // if hatches intersect the x axis, sort in order of ascending y coordinate to line them up in mark/jump order
				}
				else {
					xAsc(tmp_isList);  // if hatches intersect the y axis, sort in order of ascending x coordinate
				}
			}
			else
			{
				dirHatch = 0;
				if (primaryHatchDir == "x") {
					yDsc(tmp_isList);  // if hatches intersect the x axis, sort in order of descending y coordinate
				}
				else {
					xDsc(tmp_isList);  // if hatches intersect the y axis, sort in order of descending x coordinate
				}
			}

			// remove any duplicate vertices
			tmp_isListNoDuplicates = eliminateDuplicateVertices(tmp_isList);  // list must be sorted in x and y (or reverse) before calling this function

			// If everything went well, tmp_isListNoDuplicates should contain an even number of vertices.
			// however, due to the combination of "close enough = intersection" and "remove duplicates within proximity range" we
			// may end up with an odd number, which won't work.  If so, use the original list (including duplicates) - but only if
			// tmp_isList itself has an even number of vertices

			if ((tmp_isListNoDuplicates.size() % 2) == 0) {
				// tmp_isListNoDuplicates contains even number of vertices.  Assume it is valid
				// and add it to master intersection vector
				isList.insert(isList.end(), tmp_isListNoDuplicates.begin(), tmp_isListNoDuplicates.end());
			}
			else {
				// tmp_isListNoDuplicates contains odd number of vertices and cannot be used.
				// see if tmp_isList is valid
				if ((tmp_isList.size() % 2) == 0) {
					// use tmp_isList (with duplicates) instead
					isList.insert(isList.end(), tmp_isList.begin(), tmp_isList.end());
				}
				// if both tmp_isListNoDuplicates and tmp_isList contain odd numbers of elements, neither can be used.  Skip this hatch line
			}
		}
		// done with this hatch line
		tmp_isList.clear();
		tmp_isListNoDuplicates.clear();

		// check if we are finished (reached a_end)
		a_curr += hSpace;
		if (hSpace > 0) {
			if (a_curr >= a_end) { finished = true; }
		}
		else {
			if (a_curr <= a_end) { finished = true; }
		}
	}

	// determine whether we have anything to output (any actual intersections)
	if (isList.size() == 0) {
		// no result!,
		P.type = "";  // indicates null result
		P.vecSg.clear();
		return P;
	}

	if (useIslandScanning)
	{
		// Seed random function by something that changes every run, ensuring we get different orders each time we run 
		// (also affects rng when it comes to point fuzzing)
		srand(time(NULL));

		// Takes isList, superimposes the grid onto it, then slices all the vectors by the grid lines.
		vector<vertex> splitVectors;
		splitByIslands(&splitVectors, &isList, &boundingBox, islandSize, islandOffset, islandOverlap, hatchAngle, layerNum);
		isList = splitVectors;

		// For debug usage; offsets all the vertices by a super small, but noticeable random amount so it's easy to check for duplicate segments on top of each other
		// for (int i = 0; i < isList.size(); i++) {
		 	// fuzzPoint(&(isList[i]));
		// }
	}

	//*******************************
	// Build segments by assigning a segmentStyle to each successive pair of vertices.
	// This process is intimately intertwined with the manner in which hatch intersections are determined (above)
	// The code below alternates between mark and jump segments, and does not support zig-zag (connected) marks
	/* NOTE: Modified to account for island scanning cases where, due to half the islands being vertical and half horizontal, a given call can result 
		in NO output vectors; this is more frequent with large island sizes. The initial code setup using an iterator always executed once due to some
		weird iterator stuff or preincrementing or something like that, I just switched to a normal for loop. */
	segment sg;
	int jump = 0;
	vSg.reserve((int)isList.size());
	for (int i = 0; i < isList.size(); i++) {
		sg.start = isList[i];
		sg.end = isList[i + 1];
		if (jump == 0)
		{
			sg.idSegStyl = hatchSegStyle;
			sg.isMark = 1;
			jump = 1;  // sets up next segment to be a jump
		}
		else
		{
			sg.idSegStyl = jumpSegStyle;
			sg.isMark = 0;
			jump = 0;
		}
		//push segment to path
		vSg.push_back(sg);
	}

	// finalize the output path
	P.vecSg = vSg;
	P.tag = rProfile.Tag;
	P.SkyWritingMode = rProfile.hatchSkywriting;
	P.type = "hatch";

	return P;
}

path hatchOPT(layer &L, vector<int> regionIndex, regionProfile &rProfile, double offset, double hatchAngle, double a_min, double a_max, bool outputIntegerIDs, vector<vertex> boundingBox)
{
	/*	L: pointer to the layer structure
	regionIndex: list of regions to be hatched.  Must all have the same tag identified in rProfile
		Should indicate all related regions *and* holes, because all outer edges will be assumed to contain surface area and not holes
	rProfile: pointer to the region profile for these regions
	offset: hatch offset (from rProfile) with all contour offsets added
	hatchAngle: hatch angle in degrees.  0 = horizontal to the left
	a_min, a_max: min/max intercepts of hatch lines with vertices in L.vList, based on whether hatchAngle is closer to x or y axis
		if hatchAngle is [315 to 45] or [135 to 225], these should be min/max y coordinates
		however if hatchAngle is [45 to 135] or [225 to 315], these should be min/max x coordinates
	outputIntegerIDs: whether to use auto-generated integer ID's for segments styles, or the original string ID's
	boundingBox: vector of min/max x and y coordinates of this layer
	*/

	path P;  // path of hatch segments to be created
	vector<edge> edgeList;

	double hSpace = rProfile.resHatch;  // hSpace is the spacing between hatches, which will be modified by hatchAngle
	double a_curr;  // starting intercept coordinate of the current hatch line

	int dirHatch = 0;  // controls hatch direction; 0 is positive hatch direction, 1 = negative hatch direction
					   // current code flips dirHatch on every hatch line so that the hatches zig-zag.
					   // OPTIONAL:  allow flipping of dirHatch to be controlled by user.  If disabled, this would cause hatches to all be drawn in one direction
					   // OPTIONAL:  in addition, all user to set the value of dirHatch.  Not so useful, since dirHatch==1 would change behavior with hatch angle

	// Determine primary hatching orientation (toward x or y axis), and adjust hSpace to account for hatchAngle
	double hatchAngle_rads = hatchAngle * 3.14159265358979323846 / 180;  // hatch angle in radians
	double hatchFunctionValue;  // sine or tangent value of hatchAngle_rads to be computed once and passed to findIntersection
	string primaryHatchDir;     // either y (progress along y axis) or x (progress along x axis), depending on hatchAngle
	// Determine whether to generate hatches progressing along the vertical or horizontal axis based on hatch angle
	if ((((int)hatchAngle + 315) % 180) > 90) {
		// hatchAngle is [315 to 45] or [135 to 225] degrees, i.e. closer to the x than y axis (and potentially parallel to the x axis)
		// we will increment hatches along the y axis, which will never be parallel to the hatch direction
		primaryHatchDir = "y";
		hSpace = hSpace / cos(hatchAngle_rads);  // this is the hatch-to-hatch spacing upon intersection with y axis
		hatchFunctionValue = tan(hatchAngle_rads);
	}
	else {
		// hatchAngle is [45 to 135] or [225 to 315] degrees, i.e. closer to the y than x axis (and potentially parallel to the y axis)
		// we will increment hatches along the x axis, which will never be parallel to the hatch direction
		primaryHatchDir = "x";
		hSpace = hSpace / sin(hatchAngle_rads);  // this is the hatch-to-hatch spacing upon intersection with x axis
		hatchFunctionValue = cos(hatchAngle_rads) / sin(hatchAngle_rads);  // can't use 1/tan, which fails at 90 degrees;
	}

	int numRegions = regionIndex.size();
	int naiveTmpSize = numRegions * 3;  // estimated number of intersections per hatch line, for sizing temporary vectors

	// Identify the segment styles, in terms of either their string or auto-generated integer ID's
	string hatchSegStyle, jumpSegStyle;
	if (outputIntegerIDs == true) {
		hatchSegStyle = to_string(rProfile.hatchStyleIntID);
		jumpSegStyle = to_string(rProfile.jumpStyleIntID);
	}
	else {
		hatchSegStyle = rProfile.hatchStyleID;
		jumpSegStyle = rProfile.jumpStyleID;
	}

	//*** Offset the edges in each region (per hatch offset) and combine into a list of offsetted edges across all regions using this tag
	vector<vector<edge>> unusedVar;  // alternate output for edgeOffset, which must exist even if it won't be populated
	edgeOffset(L, regionIndex, edgeList, unusedVar, offset, false);  // false=return all offset edges in one (jumbled) vector
	
	// determine whether we have anything to output (any polygons which survived offsetting)
	if (edgeList.size() == 0) {
		// no result!
		P.type = "";  // indicates null result
		P.vecSg.clear();
		return P;
	}

	// *** Identify intersection points between edges and parallel-spaced hatch lines by stepping from a_min to a_max in increments of hSpace
	double a_start, a_end;
	if (hSpace > 0) {
		a_start = a_min;
		a_end = a_max;
	}
	else {
		a_start = a_max;
		a_end = a_min;
	}
	bool finished = false;  // use this to determine when we've reached a_end, since we could be incrementing or decrementing via hSpace
	// We start at a_start + hSpace because a hatch exactly at a_start would have no length (it would be a one-point hatch)
	a_curr = a_start + hSpace;  // a_curr is the x or y intersection of the current hatch line

	vector<hRegion> hRegionList;
	vector<hRegion> tmp_hRegionList;
	int isNum_curr, isNum_prev;
	int hStart = 0;

	//get intersection points between edges and hatch lines
	while (!finished)
	{
		vector<vertex> tmp_isList, tmp_isListNoDuplicates;  // temporary list of intersections between one hatch line and all edges for this region tag
		for (vector<edge>::iterator it = edgeList.begin(); it != edgeList.end(); ++it)
		{
			//get intersection with edge
			vertex is;
			int isB;
			isB = findIntersection(&is, hatchAngle, hatchAngle_rads, boundingBox, a_curr, *it, hatchFunctionValue);
			
			// if an intersection was found (return value == 0), add to tmp intersection vector
			if (isB == 0)
			{
				tmp_isList.push_back(is);
			}	
		}

		// Now that we've evaluated a particular hatch line against all edges, process whatever intersections we have
		//
		// If there are intersections, sort them, then eliminate any duplicates
		if (tmp_isList.size() > 0) {
			// Sort tmp_isList in x or y-direction, according to hatch direction
			//
			// Note - this method will always begin with intersection with an outer edge and proceed across
			// an area to be marked.  No matter how many intersections are crossed, the tmp_isList
			// will alternate between mark and jump segments.
			// For instance, given six vertices in tmp_isList (0-5), the resulting edges will form
			// segments 0-1 (mark), 1-2 (jump), 2-3 (mark), 3-4 (jump), 4-5 (mark), followed by a jump
			// to the next vertex, if any
			if (dirHatch == 0)
			{
				dirHatch = 1;
				if (primaryHatchDir == "x") {
					yAsc(tmp_isList);  // if hatches intersect the x axis, sort in order of ascending y coordinate to line them up in mark/jump order
				}
				else {
					xAsc(tmp_isList);  // if hatches intersect the y axis, sort in order of ascending x coordinate
				}
			}
			else
			{
				dirHatch = 0;
				if (primaryHatchDir == "x") {
					yDsc(tmp_isList);  // if hatches intersect the x axis, sort in order of descending y coordinate
				}
				else {
					xDsc(tmp_isList);  // if hatches intersect the y axis, sort in order of descending x coordinate
				}
			}

			// remove any duplicate vertices
			tmp_isListNoDuplicates = eliminateDuplicateVertices(tmp_isList);  // list must be sorted in x and y (or reverse) before calling this function

			// If everything went well, tmp_isListNoDuplicates should contain an even number of vertices.
			// however, due to the combination of "close enough = intersection" and "remove duplicates within proximity range" we
			// may end up with an odd number, which won't work.  If so, use the original list (including duplicates) - but only if
			// tmp_isList itself has an even number of vertices
			if ((tmp_isListNoDuplicates.size() % 2) == 0) {
				// tmp_isListNoDuplicates contains even number of vertices.  Assume it is valid and use it as the intersection list
				tmp_isList = tmp_isListNoDuplicates;
			}
			else {
				// tmp_isListNoDuplicates contains odd number of vertices and cannot be used.
				// see if tmp_isList is valid instead
				if ((tmp_isList.size() % 2) != 0) {
					// tmp_isList (with duplicates) is also not valid.  Can't use either of these, so eliminate all intersections
					tmp_isList.clear();
				}
			}
		}

		// at this point tmp_isList contains an sorted list of intersections
		hRegion hrg;
		segment sg;
		//create subregions to divide and conquer where each subregion has no internal voids

		//first hatch
		if (hStart == 0)
		{
			isNum_curr = tmp_isList.size();
			for (int i = 0; i < isNum_curr; i = i + 2)
			{
				hrg.start = tmp_isList[i];
				hrg.end = tmp_isList[i + 1];
				sg.start = tmp_isList[i];
				sg.end = tmp_isList[i + 1];
				sg.idSegStyl = hatchSegStyle;
				sg.isMark = 1;
				(hrg.vecSg).push_back(sg);
				tmp_hRegionList.push_back(hrg);
			}
			hStart = 1;
		}
		else
		{
			isNum_prev = isNum_curr;
			isNum_curr = tmp_isList.size();

			//if no new regions
			if (isNum_curr == isNum_prev)
			{
				int i = 0;
				for (vector<hRegion>::iterator it = tmp_hRegionList.begin(); it != tmp_hRegionList.end(); ++it)
				{
					vertex tmp_v;
					tmp_v = (*it).end;
					sg.start = (tmp_v);
					if (dirHatch)
						sg.end = (tmp_isList[i ]);
					else
						sg.end = (tmp_isList[i]);
					sg.idSegStyl = jumpSegStyle;
					sg.isMark = 0;
					((*it).vecSg).push_back(sg);
					if (dirHatch)
					{
						sg.start = (tmp_isList[i +1]);
						sg.end = (tmp_isList[i]);
						(*it).end = (tmp_isList[i]);
					}
					else
					{
						sg.start = (tmp_isList[i]);
						sg.end = (tmp_isList[i + 1]);
						(*it).end = (tmp_isList[i + 1]);
					}
					sg.idSegStyl = hatchSegStyle;
					sg.isMark = 1;
					((*it).vecSg).push_back(sg);
					i = i + 2;
				}
			}

			// if new regions encountered
			else
			{	
				hRegionList.insert(hRegionList.end(), tmp_hRegionList.begin(), tmp_hRegionList.end());
				vector<hRegion> tmp_hRegionList_old;
				tmp_hRegionList_old = tmp_hRegionList;
				tmp_hRegionList.clear();
				tmp_hRegionList.shrink_to_fit();
				for (int i = 0; i < isNum_curr; i = i + 2)
				{
					if (dirHatch)
					{
						hrg.start = (tmp_isList[i]);
						hrg.end = (tmp_isList[i+1]);
						sg.start = (tmp_isList[i]);
						sg.end = (tmp_isList[i+1]);
					}
					else
					{
						hrg.start = (tmp_isList[i]);
						hrg.end = (tmp_isList[i+1]);
						sg.start = (tmp_isList[i]);
						sg.end = (tmp_isList[i+1]);
					}
					sg.idSegStyl = hatchSegStyle;
					sg.isMark = 1;
					(hrg.vecSg).push_back(sg);
					tmp_hRegionList.push_back(hrg);
				}
			}
		}
		// check if we are finished (reached a_end)
		a_curr += hSpace;
		if (hSpace > 0) {
			if (a_curr >= a_end) { finished = true; }
		}
		else {
			if (a_curr <= a_end) { finished = true; }
		}
	}
	
	// determine whether we have anything to output (any actual intersections)
	if (hRegionList.size() == 0) {
		// no result!
		P.type = "";  // indicates null result
		P.vecSg.clear();
		return P;
	}

	//all hatches have now been generated.  put all disjoints regions in one vector
	hRegionList.insert(hRegionList.end(), tmp_hRegionList.begin(), tmp_hRegionList.end());
	vector<segment> vsg;
	vector<int> mapHrg(hRegionList.size(), 0);
	int cnt = hRegionList.size();
	vector<int> optHrg;
	int curr = 0;

	//nearest neighbor algorithm for connecting the regions
	while (cnt>0)
	{
		optHrg.push_back(curr);
		mapHrg[curr] = 1;
		cnt--;
		vertex v1 = (hRegionList[curr]).end;
		int tmp_opt = 0;
		double cost = (double)MAXCOST;
		for (int i = 0; i < (int)mapHrg.size(); i++)
		{
			if (mapHrg[i] == 0)
			{
				vertex v2 = (hRegionList[i]).start;
				double tmp_cost = dist(v1, v2);
				if (tmp_cost < cost)
				{
					cost = tmp_cost;
					tmp_opt = i;
				}
			}
		}
		curr = tmp_opt;
	}
	hRegion hrg_s, hrg_f;
	if (optHrg.size()>1)
	{
		for (vector<int>::iterator it = optHrg.begin(); it != optHrg.end() - 1; ++it)
		{
			segment sg;
			hrg_s = hRegionList[*it];
			hrg_f = hRegionList[*(it + 1)];
			sg.start = hrg_s.end;
			sg.end = hrg_f.start;
			sg.idSegStyl = jumpSegStyle;
			sg.isMark = 0;
			(hrg_s.vecSg).push_back(sg);
			vsg.insert(vsg.end(), (hrg_s.vecSg).begin(), (hrg_s.vecSg).end());
		}
	}
	else
	{
		hrg_f = hRegionList[optHrg[0]];
	}
	vsg.insert(vsg.end(), (hrg_f.vecSg).begin(), (hrg_f.vecSg).end());

	////////////////
	//Patch to make sure there are no gaps in the laser path; should not be required if the algorithm works correctly
	vsg.insert(vsg.end(), (hrg_f.vecSg).begin(), (hrg_f.vecSg).end());
	vertex vEnd, vStart, vn, tmpv;
	vector<segment> segVecNoHoles;
	segment sg;
	for (vector<segment>::iterator iSeg = vsg.begin(); iSeg != vsg.end() - 1; ++iSeg)
	{
		segVecNoHoles.push_back(*iSeg);
		vEnd = (*iSeg).end;
		vStart = (*(iSeg + 1)).start;
		vn = (*(iSeg + 1)).end;
		if ((vEnd.x != vStart.x) || (vEnd.y != vStart.y))
		{
			if (vEnd.x == vn.x && vEnd.y == vn.y)
			{
				tmpv = (*(iSeg + 1)).end;
				(*(iSeg + 1)).end = (*(iSeg + 1)).start;
				(*(iSeg + 1)).start = tmpv;
			}
			else 
			{
				sg.start = vEnd;
				sg.end = vStart;
				sg.idSegStyl = jumpSegStyle;
				sg.isMark = 0;
				segVecNoHoles.push_back(sg);
			}
		}
	}
	
	// finalize the output path
	P.vecSg = segVecNoHoles;
	P.tag = rProfile.Tag;
	P.SkyWritingMode = rProfile.hatchSkywriting;
	P.type = "hatch";

	return P;
}

path contour(layer &L, vector<int> regionIndex, regionProfile &rProfile, double offset, vector<vertex> BB, bool outputIntegerIDs)
{
	/*	L is a pointer to the layer, which contains regions, edges and vertices.
		regionIndex is a vector of region numbers (as measured from the start of (*L).s.rList) to be contoured together
		regionProfile* is a pointer to the region profile to be used
		offset is the current contour offset, which increases with the number of contours drawn.  Positive offset pulls edges "inward" (negative offsets not allowed)
		BB is the bounding box of the layer L, used to determine if any point lands out of bounds

		Process:
			Iterate over regions in regionProfile
				Offset the individual region
				If this is not the first region, add a jump from prior region endpoint to first coordinate of this region
				Iterate over offset edges in this region.  Add each edge to segment vector
				Save the final segment endpoint into a temporary variable
			Create a path out of the accumulated segments
	*/

	// Identify the values we'll need
	path P;
	string tag = rProfile.Tag;
	string markSegStyle, jumpSegStyle;
	if (outputIntegerIDs == true) {
		markSegStyle = to_string(rProfile.contourStyleIntID);
		jumpSegStyle = to_string(rProfile.jumpStyleIntID);
	}
	else {
		markSegStyle = rProfile.contourStyleID;
		jumpSegStyle = rProfile.jumpStyleID;
	}

	vector<segment> vSg;  // holds the output contour segments
	vector<vector<edge>> allOffsetEdges;  // holds return from edgeOffset containing individual vectors for each offsetted region
	segment sg;
	segment jumpConnector;  // holds endpoint of prior region so that a jump can be added between regions
	jumpConnector.start.x = 0.0;  // initialize jumpConnector
	jumpConnector.start.y = 0.0;
	jumpConnector.end = jumpConnector.start;
	jumpConnector.idSegStyl = jumpSegStyle;
	jumpConnector.isMark = 0;

	// Estimate the number of segments to enable pre-sizing the output
	int numSegments = 0;
	for (int r = 0; r != regionIndex.size(); ++r)
	{
		numSegments += L.s.rList[regionIndex[r]].eList.size() + 2;  // +2 for a segment to close the contour and jump to next region, if any
	}
	vSg.reserve(numSegments);

	//*** Offset all regions (under this tag) at once, and get a vector of individual region vectors
	vector<edge> tmpEdgeList;  // will not be used, but need to define such a variable for edgeOffset
	edgeOffset(L, regionIndex, tmpEdgeList, allOffsetEdges, offset, true);  // true=return output in allOffsetEdges, containing a vector for each region

	// determine whether we have anything to output (any polygons which survived offsetting)
	if (allOffsetEdges.size() == 0) {
		// no result!
		P.type = "";  // indicates null result
		P.vecSg.clear();
		return P;
	}

	// Iterate over regions in allOffsetEdges
	int numPopulatedRegions = 0;  // counts the number of regions with non-zero segments after offsetting
	for (int r = 0; r != allOffsetEdges.size(); ++r)
	{
		if (allOffsetEdges[r].size() > 0) {
			// only process this offsetted region if there is something there
			numPopulatedRegions += 1;

			// If this is not the first region we're processing, add a jump from prior region endpoint to first coordinate of this region
			if (numPopulatedRegions > 1)
			{
				jumpConnector.end = (*(allOffsetEdges[r]).begin()).s;
				vSg.push_back(jumpConnector);
			}

			// iterate over offsetEdgeList and create contour segments
			for (vector<edge>::iterator e = (allOffsetEdges[r]).begin(); e != (allOffsetEdges[r]).end(); ++e)
			{
				sg.start = (*e).s;
				sg.end = (*e).f;
				sg.idSegStyl = markSegStyle;
				sg.isMark = 1;
				if (findInt(BB, (*e).s, (*e).f) > 0) // check failsafe - if vertex is out of bounds, replace with a jump
				{
					sg.idSegStyl = jumpSegStyle;
					sg.isMark = 0;
				}
				vSg.push_back(sg);
			}

			// Save the endpoint from this region into a temporary variable
			jumpConnector.start = sg.end;
		}
	}

	// Convert the segments into a path structure
	P.vecSg = vSg;  // will have size zero if no result
	P.tag = tag;
	P.type = "contour";
	P.SkyWritingMode = rProfile.cntrSkywriting;
	if (vSg.size() == 0) {
		P.type = "";  // indicates no result
		P.vecSg.clear();
	}

	return P;
}
