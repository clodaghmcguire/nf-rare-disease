#!/usr/bin/awk -f
BEGIN { OFS="\t" }
{
	if ($1 ~ /(browser|track)/)
	{ }
	else if (NF<4) 
	{
		print $1,$2,$3,"SLICE_"NR 
	}
	else if (NF>=4)
	{ 
		print $1,$2,$3,$4
	}
}
