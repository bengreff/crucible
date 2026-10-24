#!/bin/sh
# Summarises the trace's clips: count, largest ratio (|rho k| over the level's largest rho k), time span, levels, cells.
O=/tmp/inq2c/mac_rfx/out
echo "trace_exit $(cat $O/trace/trace_exit)"
F=$O/trace/clips.txt
grep -c "trace clip" $F
grep "trace clip" $F | awk '{t=$4; nz=$7; r=$NF+0; n++; if(r>mx){mx=r; line=$0}; if(n==1) t0=t; t1=t; lv[nz]++; c[nz" "$9" "$11]++}
  END{print "clips",n,"t from",t0,"to",t1; print "largest ratio",mx; print line; for(k in lv) print "onto nz",k,lv[k]; m=0; for(k in c) m++; print "distinct cells",m}'
grep "trace clip" $F | awk '{print $NF}' | sort -g | awk '{a[NR]=$1} END{print "median ratio",a[int(NR/2)+1],"count above 1e-6:"; }'
grep "trace clip" $F | awk '$NF+0>1e-6' | wc -l
grep "trace clip" $F | awk '{print $9, $11}' | sort | uniq -c | sort -rn | head -5
grep -v "trace clip" $F | grep -E "maximum resident|real" | head -3
grep -E "could not|error" $O/trace/tr.log | head -3
