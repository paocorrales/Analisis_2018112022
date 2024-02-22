#!/bin/bash

module load R

FCST='2018112200'	# Fecha de pronóstico con formato YYYYMMDDHH
EXP='E8'		# Nombre del experimento
TMEM=60                 # Cantidad de miembros del ensamble

MEM=1

while [ $MEM -le $TMEM ]
do
	echo "************** Member #"${MEM}" ***************" 

	Rscript interp_sondeos.R $EXP $(printf '%02d' $((10#$MEM))) $FCST 
	MEM=$[$MEM+1]
done

