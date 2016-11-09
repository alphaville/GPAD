function  fid = gen_print_mat(fid, A, name, size_variable, is_const)
sizeA = numel(A);
if (nargin<=4 || (nargin==5 && isempty(is_const))), type = 'creal_t'; else type = is_const; end
if (nargin<=3 || (nargin>=4 && isempty(size_variable)))
    fprintf(fid, '%s ', type);
    fprintf(fid, name);
    fprintf(fid, '[%d]={', sizeA);
else
    fprintf(fid, type);
    fprintf(fid, name);
    fprintf(fid, '[%d]={', size_variable);
end
str_len = 0;
if (sum(sum(A~=0))==0)
    fprintf(fid, '%g', 0.0);
else
for i=1:sizeA
    str =num2str(A(i));
    fprintf(fid, '%g', A(i));
    str_len = str_len + numel(str) + 2;
    if (i~=sizeA), fprintf(fid, ', '); end
    if (mod(i,100)==0), fprintf(fid, '\n\t\t'); str_len = -8; end
end
end
fprintf(fid, '};\n');
end