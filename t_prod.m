function a = t_prod(b, c)
    % tubal product of 2 3d-tensors
    if ndims(b) ~= 3 || ndims(c) ~= 3
       error("dimension of tensors should be equal 3")
    end
    
    [n1, n2, n3] = size(b);
    [tmp2, n4, tmp3] = size(c);
    
    if n2 ~= tmp2 || n3 ~= tmp3
        error("Not suitable dimensions")
    end
    
    a = zeros(n1, n4, n3);
    b = fft(b, [], 3);
    c = fft(c, [], 3);
    halfn3 = round(n3/2);
    a(:,:,1) = b(:,:,1) * c(:,:,1);
    for i=2:halfn3
        a(:,:,i) = b(:,:,i) * c(:,:,i);
        a(:,:,n3+2-i) = conj(a(:,:,i));
    end
    
    if mod(n3,2) == 0
        i = halfn3+1;
        a(:,:,i) = b(:,:,i)*c(:,:,i);
    end
    a = ifft(a, [], 3);

end

