function I_t = squarePulse(t,Imax,dur)

if(0 <= t < dur)
    I_t = Imax;
else
    I_t = 0;
end