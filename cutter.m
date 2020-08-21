%% Define helper funcion cutter:
function eps = cutter(self,locs,shift_left,lstim,temp_len)
    eps = zeros(lstim-1,temp_len);
      for i=1:lstim-1
        diff=locs(i+1)-locs(i);
        eps(i,1:diff)=self(locs(i)-shift_left:locs(i)+diff-shift_left-1);
      end
end