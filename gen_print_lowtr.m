function  fid = gen_print_lowtr(fid, A, name, size_variable)
nRows = size(A,1);
sizeA = (nRows*(nRows+1))/2;
if (nargin<=3)
    fprintf(fid, 'creal_t ');
    fprintf(fid, name);
    fprintf(fid, '[%d]={', sizeA);
else
    fprintf(fid, 'creal_t ');
    fprintf(fid, name);
    fprintf(fid, '[%d]={', size_variable);
end
index = 0;
    for i=1:nRows
        for j=1:i
            index = index + 1;
            fprintf(fid, '%.16f', A(i,j));  
            if (index<sizeA), fprintf(fid, ', '); end
            if (mod(index,10)==0), fprintf(fid, '\n\t\t'); str_len = -8; end
        end
    end
    fprintf(fid, '};\n');
end