function out = correlationvec(in1,in2)
%correlation between 2 vecotre 
out = (dot(in1-mean(in1),in2-mean(in2)))/(norm(in1-mean(in1))*(norm(in2-mean(in2))));
end