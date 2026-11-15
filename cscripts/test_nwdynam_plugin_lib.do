* helpers shared by the plugin-versus-reference tests (random event sequences)
mata:
mata set matastrict off
real matrix mkev(real scalar n, real scalar N, real scalar ties, real scalar valued) {
	real matrix ev
	real scalar i, a, b, t, v
	ev = J(N,4,1); t = 0
	for (i=1;i<=N;i++) {
		a = runiformint(1,1,1,n); b = runiformint(1,1,1,n); while (b==a) b = runiformint(1,1,1,n)
		if (!(ties & runiform(1,1) < 0.3)) t = t + 0.1 + runiform(1,1)*2
		v = 1
		if (valued == 1) v = 1 + runiformint(1,1,0,3)
		if (valued == 2) v = runiformint(1,1,-1,2)
		ev[i,.] = (a,b,t,v)
	}
	return(ev)
}
end
