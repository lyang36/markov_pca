function A = createRandLumpGraph(vertNum, rk, p)
   %%create a random lumpable graph:
   %%the graph contains vertNum nodes but only contains rk meta-nodes
   %%vertNum should be a multiplier of rk
   %%p -- the probability of generating an edge
   
   mat0 = rand(rk);
   mat0 = ((mat0 + mat0) / 2) < p;
   for i = 1:rk
       mat0(i) = 0;
   end
   
   mat1 = zeros(vertNum);
   b = vertNum / rk;
   
   %fullBiGraph = [zeros(b), ones(b); ones(b), zeros(b)];
    
   for i = 1:rk
       for j = 1:rk
           if(mat0(i, j) == 1)
                pr = rand();
                mat1(((i - 1) * b + 1):((i - 1) * b + b), ((j - 1) * b + 1):((j - 1) * b + b)) = ones(b) * pr;
                mat1(((j - 1) * b + 1):((j - 1) * b + b), ((i - 1) * b + 1):((i - 1) * b + b)) = ones(b) * pr;
                if(i == j)
                    mat1(((j - 1) * b + 1):((j - 1) * b + b), ((i - 1) * b + 1):((i - 1) * b + b)) = (ones(b) - eye(b)) * pr;
                end
           end
       end
   end
   degs = sum(mat1);
   %tot_degs = sum(degs);
   A = diag(degs) \ mat1;
end

   


