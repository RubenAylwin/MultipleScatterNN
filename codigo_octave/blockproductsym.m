function z = blockproductsym(A,v,N,Nrows)

	z = zeros(size(v));

  for ii=1:Nrows


		for jj=1:(ii-1)




       z((jj-1)*N+(1:N)) = z((jj-1)*N+(1:N)) +A{jj+ii*(ii-1)/2}.'*v((ii-1)*N+(1:N));
       z((ii-1)*N+(1:N)) = z((ii-1)*N+(1:N))+A{jj+ii*(ii-1)/2}*v((jj-1)*N+(1:N));


    end

    z((ii-1)*N+(1:N)) = z((ii-1)*N+(1:N))+A{ii+ii*(ii-1)/2}*v((ii-1)*N+(1:N));



  end
end


##function z = blockproductsym(A,v)
##
##  Nblocks = size(v,1);
##
##	z = cell(size(v));
##
##	for ii=1:Nblocks
##
##     z{ii} = zeros(size(v{ii}));
##
##  end
##
##  for ii=1:Nblocks
##
##
##		for jj=1:(ii-1)
##
##
##
##
##       z{jj} = z{jj} +A{jj+ii*(ii-1)/2}'*v{ii};
##       z{ii} = z{ii}+A{jj+ii*(ii-1)/2}*v{jj};
##
##
##    end
##
##    z{ii} = z{ii}+A{ii+ii*(ii-1)/2}*v{ii};
##
##
##
##  end
##end
