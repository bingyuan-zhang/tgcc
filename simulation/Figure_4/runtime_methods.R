# CARP 100/200/500/1000
# SLC 100/200/500/1000/5000/10000
# SSNAL 100/200/500/1000/5000/10000
# CPAINT 100/200/500/1000/5000/10000/500000
# TGCC 100/200/500/1000/5000/10000/500000

carpto = 4
Tcarp = data.frame(method = rep("CARP", carpto), sample = nsize[1:carpto], time = rep(0,carpto))
Tcarptotal = data.frame(method = rep("CARPtotal", carpto), sample = nsize[1:carpto], time = rep(0,carpto))
for(i in 1:carpto) {
  t = 0
  tall = 0
  for(j in 1:Repeat) {
    print(c(i,j))
    weights = sparse_rbf_kernel_weights()
    t0 = proc.time()
    w = weights(Datalist[[i]])
    tw = proc.time()[3] - t0[3]

    carp.fit = clustRviz::CARP(Datalist[[i]])
    tall = tall + carp.fit$time
    t = t + carp.fit$fit_time + tw
  }
  Tcarp$time[i] = t/Repeat
  Tcarptotal$time[i] = tall/Repeat
}

slcto = 6
Tslc = data.frame(method = rep("SLC", slcto), sample = nsize[1:slcto], time = rep(0,slcto))
for(i in 1:slcto) {
  print(i)
  t = 0
  for(j in 1:Repeat) {
    tstart = proc.time()
    res = hclust(d = dist(Datalist[[i]]), method = "single")
    tend = proc.time()
    t = t + (tend - tstart)[3]
  }
  Tslc$time[i] = t/Repeat
}

cpaintto = 10
Tcpaint = data.frame(method = rep("CPAINT", cpaintto), sample = nsize[1:cpaintto], time = rep(0,cpaintto))
for(i in 1:cpaintto) {
  gc()
  t = 0
  for(j in 1:Repeat) {
    print(c(i,j))
    tstart = proc.time()
    lam_max = find_lambda(Datalist[[i]]);
    Lam = seq(1, lam_max, length.out = 50)
    res = cpaint(Datalist[[i]],Lam)
    tend = proc.time()
    t = t + (tend - tstart)[3]
  }
  Tcpaint$time[i] = t/Repeat
}
