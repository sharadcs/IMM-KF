function G = input_switch(t,ts,A_beg,A_end,B_beg,B_end)

if exist('B_end','var')
    %    ts=0.25;
    %A_end
    %A_beg
    G = inv(A_end)*(expm(A_end*(ts-t))-eye(4))*B_end ...
      + expm(A_end*(ts-t))*inv(A_beg)*(expm(A_beg*(t))-eye(4))*B_beg;
else
    %    ts=0.0005;
    G = inv(A_end)*(expm(A_end*(ts-t))-eye(size(A_end,1)))*B_beg ...
      + expm(A_end*(ts-t))*inv(A_beg)*(expm(A_beg*(t))-eye(size(A_end,1)))*B_beg;
end

end
