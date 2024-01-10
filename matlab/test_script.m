function test_script(thedir, k)

t0=datestr(now); 
pause(5)
t1=datestr(now);
save(sprintf('%s/test_%d', thedir, k), 't0','t1');