function[min ,max ]= f_calc_minmax( min, max ,type , s , id , tp)
if nargin < 6
    tp = ':';
end
sim_string = strcat( 's.',type,'values(',tp ,', id)' );
eval(strcat('max( max < ', sim_string,')=  ', 's.', type,'values( max <', sim_string,', id );'))
eval(strcat('min( min > ', sim_string,')=  ', 's.', type,'values( min >', sim_string,', id );'))

end