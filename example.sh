make 

./tgsoa \
  --problem=etc/sphere_rand_ol_100x1.tsp

./tgsoa \
  --problem=etc/sphere_rand_ol_100x1.tsp \
  --lkh 1

./tgsoa \
  --problem=etc/sphere_rand_ol_100x1.tsp \
  --glkh 1 \
  --glkh-instance=etc/gtsp/sphere_rand_ol_100x10.gtsp 