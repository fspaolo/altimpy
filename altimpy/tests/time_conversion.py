import datetime as dt
import altimpy as ap

t1 = dt.datetime(1, 1, 1, 0, 0, 0)
t2 = dt.datetime(1, 1, 1, 0, 0, 1)
t3 = dt.datetime(1, 1, 1, 0, 1, 0)
t4 = dt.datetime(1, 1, 1, 1, 0, 0)
t5 = dt.datetime(2013, 1, 1, 0, 0, 0)
t6 = dt.datetime(2013, 1, 1, 0, 0, 1)
t7 = dt.datetime(2013, 1, 1, 0, 1, 0)
t8 = dt.datetime(2013, 1, 1, 1, 0, 0)

y1 = ap.date2year(t1)
y2 = ap.date2year(t2)
y3 = ap.date2year(t3)
y4 = ap.date2year(t4)
y5 = ap.date2year(t5)
y6 = ap.date2year(t6)
y7 = ap.date2year(t7)
y8 = ap.date2year(t8)

d1 = ap.year2date(y1)
d2 = ap.year2date(y2)
d3 = ap.year2date(y3)
d4 = ap.year2date(y4)
d5 = ap.year2date(y5)
d6 = ap.year2date(y6)
d7 = ap.year2date(y7)
d8 = ap.year2date(y8)

print '(def date -> date2year -> year2date)'
print 'original_date    converted_date      converted_year'
print t1, ' ', d1[0], ' ', y1[0]
print t2, ' ', d2[0], ' ', y2[0]
print t3, ' ', d3[0], ' ', y3[0]
print t4, ' ', d4[0], ' ', y4[0]
print t5, ' ', d5[0], ' ', y5[0]
print t6, ' ', d6[0], ' ', y6[0]
print t7, ' ', d7[0], ' ', y7[0]
print t8, ' ', d8[0], ' ', y8[0]

t1 = 1
t2 = 1.00000003171
t3 = 1.00000190259
t4 = 1.00011415525
t5 = 2013
t6 = 2013.00000003
t7 = 2013.0000019
t8 = 2013.00011416

d1 = ap.year2date(t1)
d2 = ap.year2date(t2)
d3 = ap.year2date(t3)
d4 = ap.year2date(t4)
d5 = ap.year2date(t5)
d6 = ap.year2date(t6)
d7 = ap.year2date(t7)
d8 = ap.year2date(t8)

y1 = ap.date2year(d1)
y2 = ap.date2year(d2)
y3 = ap.date2year(d3)
y4 = ap.date2year(d4)
y5 = ap.date2year(d5)
y6 = ap.date2year(d6)
y7 = ap.date2year(d7)
y8 = ap.date2year(d8)

print ''
print '(def year -> year2date -> date2year)'
print 'original_year    converted_year      converted_date'
print t1, ' ', y1[0], ' ', d1[0]
print t2, ' ', y2[0], ' ', d2[0]
print t3, ' ', y3[0], ' ', d3[0]
print t4, ' ', y4[0], ' ', d4[0]
print t5, ' ', y5[0], ' ', d5[0]
print t6, ' ', y6[0], ' ', d6[0]
print t7, ' ', y7[0], ' ', d7[0]
print t8, ' ', y8[0], ' ', d8[0]
