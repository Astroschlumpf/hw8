#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

// ************************* Hilfsfunktionen *************************** //

double doInvSq(double *x){
  return 1/(sqrt(x[0]*x[0] + x[1]*x[1]));
}

double EnVal(double *p, double *q){
  return 0.5 * (p[0]*p[0] + p[1]*p[1]) - doInvSq(q);
}

// ************************ DGL-System aus Ham. ************************* //

void KepDGL(double *p, double *q, const double &dt){
  double altP[2];
  altP[0] = p[0]; altP[1] = p[1]; // dummies
  double invD = doInvSq(q); // "1/()"-Term
  
  // cout << invD << endl; // debug
  
  p[0] = p[0] - dt * invD*invD*invD * q[0];                    //
  p[1] = p[1] - dt * invD*invD*invD * q[1];                    // ausgewertete Gleichungen fuer
  q[0] = q[0] + dt * altP[0] - dt*dt * invD*invD*invD * q[0];  // p_n+1 und q_n+1
  q[1] = q[1] + dt * altP[1] - dt*dt * invD*invD*invD * q[1];  //

  //delete[] altP;
}

// ********************************************************************* //

int main(void){
  double t = 0.0; // Start: 0.0, Ende: 20*pi
  double tEnd = 20 * M_PI;
  const double dt = 0.0005; // dt = 0.05;
  int ind = 1; // Schleifenzaehler: startet mit zweitem Eintrag, erster ist (p,q)[0]
  const int iMAX = 125666; // 20*pi/.0005 + Puffer
  
  double q[2]; double p[2];
  q[0] = 0.4; q[1] = 0.; p[0] = 0.; p[1] = 2.; // sqrt(1.6/0.4) = sqrt(4) == 2

  ofstream orb("orbit.csv");
  ofstream ener("energie.csv");
  orb  << t << "\t" << q[0] << "\t" << q[1] << endl; // Startposition fuer Orbit
  ener << t << "\t" << EnVal(p,q) << endl; // erster Energiewert
  // cout << t << "\t" << EnVal(p,q) << endl; // debug
  
  while((t < tEnd) && (ind < iMAX)){
    KepDGL(p,q,dt);
    
    orb  << t << "\t" << q[0] << "\t" << q[1] << endl; // Ausgabe der Orbitalkoordinaten
    ener << t << "\t" << EnVal(p,q) << endl; // Ausgabe der aktuellen Gesamtenergie
    // cout << t << "\t" << EnVal(p,q) << endl; // debug
    
    ind++;
    t += dt;
  }
  
  ener.close();
  orb.close();
  
  //delete[] q;
  //delete[] p;
  return 0;
}