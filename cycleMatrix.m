function a = cycleMatrix(len, z)
    m = zeros(len, len);
    for j = 1:len
       m(j, mod(j + 1, len) + 1) = -z;
       m( mod(j + 1, len) + 1, j) = -z;
       m( mod(j - 1, len) + 1, j) = -z;
       m(j, mod(j + len/2, len) + 1) = 1;
       m( mod(j + len/2, len) + 1, j) = 1;
       m(j, j) = 2 * z - 2;
    end
    a = m;
end