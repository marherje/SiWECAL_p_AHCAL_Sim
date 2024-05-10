# Patch to remove the smoothing bug START                                                                                                                                                                   
delete_test = ['Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 2','Error in <TH1F::Smooth>: Smooth only supported for histograms with >= 3 bins. Nbins = 1']
f = open('example.txt','r')
lst = []
for line in f:
    for word in delete_test:
        if word in line:
            line = line.replace(word,'')
    if line != '\n':
        lst.append(line)
f.close()
f = open('example.txt','w')
for line in lst:
    f.write(line)
f.close()
# Patch to remove the smoothing bug END
