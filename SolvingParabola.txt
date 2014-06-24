 1210 //=========================================================================
 1211 // Solve parabola using Cramer's rule 
 1212 //========================================================================
 1213 void PrSeedingXLayers::solveParabola(const PrHit* hit1, const PrHit* hit2, const PrHit* hit3, float& a, float& b, float& c){
 1214   
 1215   const float z1 = hit1->z() - m_geoTool->zReference();
 1216   const float z2 = hit2->z() - m_geoTool->zReference();
 1217   const float z3 = hit3->z() - m_geoTool->zReference();
 1218   
 1219   const float x1 = hit1->x();
 1220   const float x2 = hit2->x();
 1221   const float x3 = hit3->x();
 1222   
 1223   const float det = (z1*z1)*z2 + z1*(z3*z3) + (z2*z2)*z3 - z2*(z3*z3) - z1*(z2*z2) - z3*(z1*z1);
 1224   
 1225   if( fabs(det) < 1e-8 ){
 1226     a = 0.0;
 1227     b = 0.0;
 1228     c = 0.0;
 1229     return;
 1230   }
 1231   
 1232   const float det1 = (x1)*z2 + z1*(x3) + (x2)*z3 - z2*(x3) - z1*(x2) - z3*(x1);
 1233   const float det2 = (z1*z1)*x2 + x1*(z3*z3) + (z2*z2)*x3 - x2*(z3*z3) - x1*(z2*z2) - x3*(z1*z1);
 1234   const float det3 = (z1*z1)*z2*x3 + z1*(z3*z3)*x2 + (z2*z2)*z3*x1 - z2*(z3*z3)*x1 - z1*(z2*z2)*x3 - z3*(z1*z1)*x2;
 1235 
 1236   a = det1/det;
 1237   b = det2/det;
 1238   c = det3/det;
 1239   
 1240   
 1241 
 1242 
 1243 
 1244 
 1245 }
 1246 
