# MBS-raw
The raw but working version of new MBS calculation program. This program is for solving the light scattering problem for solid faceted particles (like ice or mineral crystals)

Examples of commands:
1.	Hexagonal column, geometrical optics, all scattering angles:
`-p 1 100 69.6 --ri 1.3116 0 –n 8 --random 300 400 --go --grid 180 45 360 --all –w 0.532 –o out`
2.	Particle from file “cube1.dat”, physical optics, backscattering direction, **absorption** enabled, the program should be closed, with scattering filter of 40°:
`--pf -p cube1.dat --rs 10 --ri 1.32 0.000019 -n 12 --random 1200 750 --po --grid 0.5 2 2 --all –w 1.064 –o out2 --abs --close --filter 40`
3.	To do multithread calculations you just need to split your calculation by orientations of beta or gamma angle. In that case you need to use "-b" or "-g" parameter. For example:
```
-p 1 100 69.6 --ri 1.3116 0 –n 8 --random 150 400 --po --grid 0.5 2 2 --all –w 1.064 –o out2 --abs --close -b 0 45 
-p 1 100 69.6 --ri 1.3116 0 –n 8 --random 150 400 --po --grid 0.5 2 2 --all –w 1.064 –o out2 --abs --close -b 45.3 90
```

More information see in [manual](https://github.com/Heart-Under-Blade/MBS-raw/blob/main/manual.pdf)
