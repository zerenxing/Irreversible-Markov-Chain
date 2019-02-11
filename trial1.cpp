double w[17][3]; // Boltzmann factors
void computeBoltzmannFactors (double J,double H,double T) {
	for (int i = -8; i <= 8; i += 4) {
		w[i + 8][0] = exp( - (i * J + 2 * H) / T);
		w[i + 8][2] = exp( - (i * J - 2 * H) / T);
	}
}
bool MetropolisStep (int sign ) {
	// choose a random spin
	int i = floor(Lx*std_rand());
	int j = floor(Ly*std_rand());
	if (s[i][j]== sign){
		// find its neighbors using periodic boundary conditions
		int iPrev = i == 0 ? Lx-1 : i-1;
		int iNext = i == Lx-1 ? 0 : i+1;
		int jPrev = j == 0 ? Ly-1 : j-1;
		int jNext = j == Ly-1 ? 0 : j+1;
		// find sum of neighbors
		int sumNeighbors = s[iPrev][j] + s[iNext][j] + s[i][jPrev] + s[i][jNext];
		int delta_ss = -2*s[i][j]*sumNeighbors;
		// ratio of Boltzmann factors
		double spin_ratio = w[delta_ss+8][1+s[i][j]];
		if (std_rand() < ratio) {
		s[i][j] = -s[i][j];
		return true;
		}
	double switch_ratio=
	if x[]
	}
	
	
}