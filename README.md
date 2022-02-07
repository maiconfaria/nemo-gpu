#adding parallel region to DO_3D
sed -i -e ' /DO_3D/i !$acc parallel' -e  '/END_3D/a !$acc end parallel ' $(find src.acc -type f | grep -v do_loop_substitute.h90 ) 
