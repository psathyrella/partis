from rpy2 import *

faithful_data = {"eruption_duration":[],
                 "waiting_time":[]}
		 
f = open('faithful.dat','r')

for row in f.readlines()[1:]: # skip the column header line
    splitrow = row[:-1].split(" ")
    faithful_data["eruption_duration"].append(float(splitrow[0]))
    faithful_data["waiting_time"].append(int(splitrow[1]))

f.close()

ed = faithful_data["eruption_duration"] 
edsummary = r.summary(ed)
print "Summary of Old Faithful eruption duration data"
for k in edsummary.keys():
    print k + ": %.3f" % edsummary[k]
print    
print "Stem-and-leaf plot of Old Faithful eruption duration data"
print r.stem(ed)

r.png('faithful_histogram.png',width=733,height=550)
r.hist(ed,r.seq(1.6, 5.2, 0.2), prob=1,col="lightgreen",
       main="Old Faithful eruptions",xlab="Eruption duration (seconds)")
r.lines(r.density(ed,bw=0.1),col="orange")
r.rug(ed)
r.dev_off()

long_ed = filter(lambda x: x > 3, ed)
r.png('faithful_ecdf.png',width=733,height=550)
r.library('stats')
r.plot(r.ecdf(long_ed), do_points=0, verticals=1, 
       main="Empirical cumulative distribution function of Old Faithful eruptions longer than 3 seconds")
x = r.seq(3,5.4,0.01)
r.lines(r.seq(3,5.4,0.01),r.pnorm(r.seq(3,5.4,0.01),mean=r.mean(long_ed), 
        sd=r.sqrt(r.var(long_ed))), lty=3, lwd=2, col="red")
r.dev_off()

r.png('faithful_qq.png',width=733,height=550)
r.par(pty="s")
r.qqnorm(long_ed,col="blue")
r.qqline(long_ed,col="red")
r.dev_off()

print
print "Shapiro-Wilks normality test of Old Faithful eruptions longer than 3 seconds"
sw = r.shapiro_test(long_ed)
print "W = %.4f" % sw['statistic']['W']
print "p-value = %.5f" % sw['p.value']

print
print "One-sample Kolmogorov-Smirnov test of Old Faithful eruptions longer than 3 seconds"
ks = r.ks_test(long_ed,"pnorm", mean=r.mean(long_ed), sd=r.sqrt(r.var(long_ed)))
print "D = %.4f" % ks['statistic']['D']
print "p-value = %.4f" % ks['p.value']
print "Alternative hypothesis: %s" % ks['alternative']
print
