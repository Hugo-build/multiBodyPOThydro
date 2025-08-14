function[p] = pressure(dat5p,count)
      p1 = dat5p(count,4)+1i*dat5p(count,5);                
      p2 = dat5p(count,6)+1i*dat5p(count,7);
      p3 = dat5p(count,8)+1i*dat5p(count,9);                
      p = [p1, p2, p3, 0, 0, 0]; 
