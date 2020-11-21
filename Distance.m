function [distance]= Distance(lnA,laA,lnB,laB)
  ra = 6378140; 
  rb = 6356755; 
  flatten = (ra -rb)/ra; 
  
  radLaA = deg2rad(laA);
  radLnA = deg2rad(lnA);
  radLaB = deg2rad(laB);
  radLnB = deg2rad(lnB);
  pa = atan(rb/ra*tan(radLaA));
  pb = atan(rb/ra*tan(radLaB));
  x = acos(sin(pa)*sin(pb)+cos(pa)*cos(pb)*cos(radLnA-radLnB));
  
  c1 = (sin(x)-x)*(sin(pa)+sin(pb))^2/cos(x/2)^2;
  c2 = (sin(x)+x)*(sin(pa)-sin(pb))^2/sin(x/2)^2;
  
  dr = flatten/8*(c1-c2);
  distance = ra*(x+dr);