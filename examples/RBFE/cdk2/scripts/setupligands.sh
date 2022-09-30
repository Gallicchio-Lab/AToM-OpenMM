for i in *.sdf ; do ( mol=${i%.sdf} ; cp ../ffgen/ffgen/${mol}/vacuum.mol2 ${mol}.mol2 ) ; done
for i in *.sdf ; do ( mol=${i%.sdf} ; cp ../ffgen/ffgen/${mol}/leaprc_header ${mol}.leaprc ) ; done
for i in *.sdf ; do ( mol=${i%.sdf} ; cp ../ffgen/ffgen/${mol}/vacuum.frcmod ${mol}.frcmod ) ; done
for i in *.sdf ; do ( mol=${i%.sdf}; python ../scripts/to_upper.py ${mol}.mol2; python ../scripts/to_upper.py ${mol}.frcmod; python ../scripts/to_upper.py ${mol}.leaprc; ); done
for i in *.sdf ; do ( mol=${i%.sdf}; cp ${mol}.mol2 ${mol}-p.mol2 ; cp ${mol}.frcmod ${mol}-p.frcmod ; cp ${mol}.leaprc  ${mol}-p.leaprc  ); done
