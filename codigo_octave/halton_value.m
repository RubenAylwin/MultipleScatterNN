function r = halton_value ( i, m )

%*****************************************************************************80
%
%% halton_value() computes an element of a Halton sequence.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    10 August 2016
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    John Halton,
%    On the efficiency of certain quasi-random sequences of points
%    in evaluating multi-dimensional integrals,
%    Numerische Mathematik,
%    Volume 2, pages 84-90, 1960.
%
%  Input:
%
%    integer I, the index of the element of the sequence.
%    0 <= I.
%
%    integer M, the spatial dimension.
%    1 <= M <= 100.
%
%  Output:
%
%    real R(M), the element of the sequence with index I.
%
## prime = [ ...
##        2;    3;    5;    7;   11;   13;   17;   19;   23;   29; ...
##       31;   37;   41;   43;   47;   53;   59;   61;   67;   71; ...
##       73;   79;   83;   89;   97;  101;  103;  107;  109;  113; ...
##      127;  131;  137;  139;  149;  151;  157;  163;  167;  173; ...
##      179;  181;  191;  193;  197;  199;  211;  223;  227;  229; ...
##      233;  239;  241;  251;  257;  263;  269;  271;  277;  281; ...
##      283;  293;  307;  311;  313;  317;  331;  337;  347;  349; ...
##      353;  359;  367;  373;  379;  383;  389;  397;  401;  409; ...
##      419;  421;  431;  433;  439;  443;  449;  457;  461;  463; ...
##      467;  479;  487;  491;  499;  503;  509;  521;  523;  541 ];

  prime = primes(10001).';

  i = abs ( floor ( i ) );

  t(1:m,1) = i;
%
%  Carry out the computation.
%
  prime_inv = zeros ( m, 1 );
  for i1 = 1 : m
    prime_inv(i1,1) = 1.0 / prime(i1,1);
  end

  r = zeros ( m, 1 );

  while ( any ( t ~= 0 ) )
    for j = 1 : m
      d = mod ( t(j), prime(j) );
      r(j) = r(j) + d * prime_inv(j);
      prime_inv(j) = prime_inv(j) / prime(j);
      t(j) = floor ( t(j) / prime(j) );
    end
  end

  return
end

