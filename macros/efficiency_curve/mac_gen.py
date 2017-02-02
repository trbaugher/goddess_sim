import string
emin = 100
emax = 2000
nevt = 100000
estep = 100
eunit = "keV"
e = emin
while e <= emax:
  print str(e)+str(eunit)+".mac"
  filename = str(e)+str(eunit)+".mac"
  file = open(filename, 'w')
  file.write("/gps/particle gamma\n")
  file.write("/gps/pos/centre 0. 0. 0. cm\n")
  file.write("/gps/ang/type iso\n")
  file.write("/gps/energy "+str(e)+" "+str(eunit)+"\n")
  file.write("/gps/List\n")
  file.write("/run/FileName out/efficiency_curve/"+str(e)+str(eunit)+".out"+"\n")
  file.write("/control/execute macros/efficiency_curve/beamon.mac\n")
  file.close()
  e = e + estep
