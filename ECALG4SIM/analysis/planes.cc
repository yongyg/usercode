
//plane : ax + by + cz + d = 0;  



// compute the plane equation from a set of 3 points.
// returns false if any of the points are co-incident and do not form a plane.
// A = point 1
// B = point 2
// C = point 3
// plane = destination for plane equation A,B,C,D
bool computePlane(const double A[],const double B[],const double C[],double plane[])
{
  bool ret = false;

  double vx = (B[0] - C[0]);
  double vy = (B[1] - C[1]);
  double vz = (B[2] - C[2]);

  double wx = (A[0] - B[0]);
  double wy = (A[1] - B[1]);
  double wz = (A[2] - B[2]);

  double vw_x = vy * wz - vz * wy;
  double vw_y = vz * wx - vx * wz;
  double vw_z = vx * wy - vy * wx;

  double mag = sqrtf((vw_x * vw_x) + (vw_y * vw_y) + (vw_z * vw_z));

  if ( mag > 0 )
    {

      mag = 1.0/mag; // compute the recipricol distance

      ret = true;

      plane[0] = vw_x * mag;
      plane[1] = vw_y * mag;
      plane[2] = vw_z * mag;
      plane[3] = 0.0 - ((plane[0]*A[0])+(plane[1]*A[1])+(plane[2]*A[2]));

    }
  return ret;
}


// inertesect a line semgent with a plane, return false if they don't intersect.
// otherwise computes and returns the intesection point 'split'
// p1 = 3d point of the start of the line semgent.
// p2 = 3d point of the end of the line segment.
// split = address to store the intersection location x,y,z.
// plane = the plane equation as four doubles A,B,C,D.

bool intersectLinePlane(const double p1[],const double p2[],double split[],const double plane[])
{

  double dp1 = p1[0]*plane[0] + p1[1]*plane[1] + p1[2]*plane[2] + plane[3];
  double dp2 = p2[0]*plane[0] + p2[1]*plane[1] + p2[2]*plane[2] + plane[3];

  if ( dp1 > 0 && dp2 > 0 ) return false;
  if ( dp1 < 0 && dp2 < 0 ) return false;

  double dir[3];

  dir[0] = p2[0] - p1[0];
  dir[1] = p2[1] - p1[1];
  dir[2] = p2[2] - p1[2];

  double dot1 = dir[0]*plane[0] + dir[1]*plane[1] + dir[2]*plane[2];
  double dot2 = dp1 - plane[3];

  double t = -(plane[3] + dot2 ) / dot1;

  split[0] = (dir[0]*t)+p1[0];
  split[1] = (dir[1]*t)+p1[1];
  split[2] = (dir[2]*t)+p1[2];

  return true;
}


bool intersectLinePlane_tolerance(const double p1[],const double p2[],double split[],const double plane[], double epsilon = 1E-4)
{

  double dp1 = p1[0]*plane[0] + p1[1]*plane[1] + p1[2]*plane[2] + plane[3];
  double dp2 = p2[0]*plane[0] + p2[1]*plane[1] + p2[2]*plane[2] + plane[3];
  
  if( fabs(dp1)> epsilon && fabs(dp2) > epsilon){ //within the tolereance, that means point is on the plane 
    
    if ( dp1 > 0 && dp2 > 0 ) return false;
    if ( dp1 < 0 && dp2 < 0 ) return false;

  }


  double dir[3];

  dir[0] = p2[0] - p1[0];
  dir[1] = p2[1] - p1[1];
  dir[2] = p2[2] - p1[2];

  double dot1 = dir[0]*plane[0] + dir[1]*plane[1] + dir[2]*plane[2];
  double dot2 = dp1 - plane[3];

  double t = -(plane[3] + dot2 ) / dot1;

  split[0] = (dir[0]*t)+p1[0];
  split[1] = (dir[1]*t)+p1[1];
  split[2] = (dir[2]*t)+p1[2];

  return true;
}



// compute the distance between a 3d point and a plane
double distToPlaneOLD(const double p[],const double plane[])
{
  return p[0]*plane[0] + p[1]*plane[1] + p[2]*plane[2] + plane[3];
}



// compute the distance between a 3d point and a plane
double distToPlane(const double p[],const double plane[])
{
  double dist = p[0]*plane[0] + p[1]*plane[1] + p[2]*plane[2] + plane[3];
  return dist / sqrt( plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2]);
}


double distTwoPoints(const double p1[],const double p2[]){
  
  return sqrt( pow(p1[0]-p2[0],2) + pow(p1[1]-p2[1],2) + pow(p1[2]-p2[2],2) );
  
}


double check4PointsInAplane(double p1[],double p2[],double p3[],double p4[]){
  double p[4];
  computePlane(p1,p2,p3,p);
  double y = p[0] * p4[0] + p[1] * p4[1] + p[2]*p4[2] + p[3];
  return y;
  
}


int pnpoly(int nvert, double vertx[], double verty[], double testx, double testy){
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
	 (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
  }
  return c;
}

// area2D_Polygon(): compute the area of a 2D polygon
//  Input:  int n = the number of vertices in the polygon
//          Point* V = an array of n+2 vertices with V[n]=V[0]
//  Return: the (float) area of the polygon
double area2D_Polygon( int n, double X[],double Y[]){
  double area = 0;
  int  i, j, k;   // indices

  if (n < 3) return 0;  // a degenerate polygon
  
  for (i=1, j=2, k=0; i<n; i++, j++, k++) {
    //double tmp = j<n? Y[j]: Y[0];
    area += X[i] * (Y[j] - Y[k]);
  }
  
  area += X[0] * (Y[1] - Y[n-1]);  // wrap-around term
  return area / 2.0;
}

double AreaofPolygon(int nvert,double X[], double Y[]){
  
  double  area=0. ;
  int     i, j=nvert-1  ;
  
  for (i=0; i<nvert; i++) {
    area+=(X[j]+X[i])*(Y[j]-Y[i]); 
    j=i; 
  }
  
  return area*.5; 
  
}

void distPreStepPointAllSides(double pos[],double pre[],int ieta, int iphi, float res[], int iz= 0, int inEB=1){
  
  if( inEB && !( ieta>=0 && ieta<=169 && iphi>=0 && iphi<=369)){
    cout<<"input ieta/iphi "<< ieta<<" "<<iphi<<endl; 
    exit(1);
  }
  
  double p[4];
  double ct[3]; 

  if(inEB==1){
    ct[0] = xctEBAll[ieta][iphi];
    ct[1] = yctEBAll[ieta][iphi];
    ct[2] = zctEBAll[ieta][iphi];
  }else{
    ct[0] = xctEEAll[iz][ieta][iphi];
    ct[1] = yctEEAll[iz][ieta][iphi];
    ct[2] = zctEEAll[iz][ieta][iphi];

    if(iz!=1 && iz!=0){
      cout<<"wrong iz "<< iz <<endl; 
      exit(1);
    }
    
  }
  
  double split[3];
  bool intersets[6]; 
  int inpolygon[6];
  
  double dist_ct[6];
  double dist_pre[6];
  double distpremin = 1E9; 
  int ind_distpremin = -1; 

  double split_new_all[6][2]; 
  for(int j=0; j<6;j++){
    for(int k=0; k<2; k++){
      split_new_all[j][k] = -9;
    }
  }
  
  for(int pl = 0; pl<=5; pl++){
    for(int j1=0; j1<4; j1++){
      
      if(inEB==1){
	p[j1] = map_planes_eb[ieta][iphi][4*pl+j1];
      }else{
	p[j1] = map_planes_ee[iz][ieta][iphi][4*pl+j1];
      }
    }
    
    double dd1 = distToPlane(pre,p);
    double dd0 = distToPlane(ct,p);
    
    if(distpremin>fabs(dd1)){
      distpremin = fabs(dd1); 
      ind_distpremin = pl;
    }

    dist_ct[pl] = dd0;
    dist_pre[pl] = dd1;
    
    //bool interset_strict = intersectLinePlane(pos,pre,split,p);
    bool interset = intersectLinePlane_tolerance(pos,pre,split,p);
    
    intersets[pl] = interset;
    inpolygon[pl] = 0; 
    
    if(interset){//check the split if inside the area
      double points[4][3];
      
      double points_new[4][2];
      points_new[0][0] =0;
      points_new[0][1] =0;
      
      for(int l=0; l<4; l++){ // 4 lines
	int ind1 = coindall[pl][l];
	
	if(inEB==1){
	  points[l][0] = coxEBAll[ieta][iphi][ind1];
	  points[l][1] = coyEBAll[ieta][iphi][ind1];
	  points[l][2] = cozEBAll[ieta][iphi][ind1];
	}else{
	  points[l][0] = coxEEAll[iz][ieta][iphi][ind1];
          points[l][1] = coyEEAll[iz][ieta][iphi][ind1];
          points[l][2] = cozEEAll[iz][ieta][iphi][ind1];
	}
	
      }
      points_new[1][0] = distTwoPoints(points[1],points[0]);
      points_new[1][1] = 0; 
      TVector3 v12(points[1][0]-points[0][0],points[1][1]-points[0][1],points[1][2]-points[0][2] );
      TVector3 v13(points[2][0]-points[0][0],points[2][1]-points[0][1],points[2][2]-points[0][2] );
      TVector3 v14(points[3][0]-points[0][0],points[3][1]-points[0][1],points[3][2]-points[0][2] );
      double theta = v12.Angle(v13);
      double d13 = distTwoPoints(points[0],points[2]);
      
      if( fabs(d13*cos(theta)-v12.Dot(v13)/v12.Mag() ) > 1E-10 ){
	cout<<"check1 "<< d13*cos(theta) <<" "<< v12.Dot(v13)/v12.Mag() <<" "<< d13*cos(theta) - v12.Dot(v13)/v12.Mag() <<endl;
	exit(1);
      }
      TVector3 v = v12.Cross(v13);
      if( v.z() ==0){
	cout<<"wrong v12Xv13 ?? "<< v.z()<<endl; 
	exit(1);
      }
      
      if( fabs( v.Mag()/v12.Mag() - d13 * sin(theta)) > 1E-10 ){
	cout<<"check2 "<< v.Mag()/v12.Mag() <<" "<< d13 * sin(theta) <<" "<<  v.Mag()/v12.Mag() - d13 * sin(theta) <<endl;
	exit(1);
      }
      
      points_new[2][0] = v12.Dot(v13)/v12.Mag();
      points_new[2][1] = v.z()/fabs(v.z()) * v.Mag()/v12.Mag(); 
      theta = v12.Angle(v14);
      v = v12.Cross(v14);
      if( v.z() ==0){
	cout<<"wrong v12Xv14 ?? "<< v.z()<<endl; 
	exit(1);
      }
      
      points_new[3][0] = v12.Dot(v14)/v12.Mag();
      points_new[3][1] = v.z()/fabs(v.z()) * v.Mag()/v12.Mag();
      
      //now check if the split
      double split_new[2];
      
      TVector3 v1s(split[0]-points[0][0],split[1]-points[0][1],split[2]-points[0][2]);
      v = v12.Cross(v1s);
      split_new[0] = v12.Dot(v1s)/v12.Mag();
      if(v.z()==0){
	split_new[1] = 0 ;
      }	else{
	split_new[1] = v.z()/fabs(v.z()) * v.Mag()/v12.Mag();
      }

      //now check split_new inside
      double vx[4];
      double vy[4];
      
      for(int n=0; n<4; n++){
	vx[n] = points_new[n][0];
	vy[n] = points_new[n][1];
      }
      int c = pnpoly(4,vx,vy,split_new[0],split_new[1]);
      inpolygon[pl] = c; 


      split_new_all[pl][0] = split_new[0];
      split_new_all[pl][1] = split_new[1];
      
      
    }///if interset
    

  }//all 6 planes
  
  //first find intersets and inside area
  vector<int> pl; 
  for(int j=0; j<6; j++){
    if( intersets[j] && inpolygon[j] ){
      pl.push_back(j);
    }
  }

  int ind; 
  if(pl.size()>1){

    double premin = 1E9; 
    ind = pl[0];
    for(int j=0;j<int(pl.size());j++){
      if(premin> fabs(dist_pre[pl[j]])){
	premin = fabs(dist_pre[pl[j]]); 
	ind = pl[j];
      }
    }
    
  }else if(pl.size()==1){
    ind = pl[0]; 
    

  }else{
    ind = ind_distpremin;
  }
  res[0] = dist_pre[ind] * dist_ct[ind]>0 ? fabs(dist_pre[ind]) : -fabs(dist_pre[ind]);  //postive means inside
  
  res[1] = ind; 
  res[2] = split_new_all[ind][0];
  res[3] = split_new_all[ind][1];
    
  
}

