function z = blockproduct(A,v,N,Nrows)

	z = zeros(size(v));

  for ii=1:Nrows


		for jj=1:Nrows




       z((ii-1)*N+(1:N)) = z((ii-1)*N+(1:N))+A{ii,jj}*v((jj-1)*N+(1:N));


    end

  end
end

