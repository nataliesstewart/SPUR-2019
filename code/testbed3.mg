generateHeckeAlgebra := function(e,n)
	PreL<v> := CyclotomicField(2*e); // Can use a different field/ring, e.g., PolynomialRing(Rationals)
	R<x> := PolynomialRing(PreL);
	q := v^2;
  
	vart := -1 - q^2 - q^3;
	L := ext<PreL | x^2-vart>;
	q := L!q;

	F<[T]> := FreeAlgebra(L,n-1);
	
	// Generate List Hecke Algebra Relations
	HRels := [];
	for i := 1 to n-1 do
	  HRels := Append(HRels,T[i]^2 - (q-1)*T[i]-q);
	  if i le n-2 then
	    HRels := Append(HRels,T[i]*T[i+1]*T[i] - T[i+1]*T[i]*T[i+1]);
	  end if;
	  for j := i+2 to n-1 do
	    HRels := Append(HRels,T[i]*T[j] - T[j]*T[i]);
	  end for;
	end for;

	// Define the Hecke Algebra A
	H := quo<F|HRels>;
	A, iso := Algebra(H); //iso is isomorphism A->H
	J := Inverse(iso);
	
	return H,A,T,q,L,J;
end function;

getReps := function(e,n,r,sp,H,A,T,q,L,J)
	vart := -1 - q^2 - q^3;
	vara := L!-1;
	varb := q;
	varc := vart*(q*vart-1);
	vard := Sqrt(vart)^3*(q+1);
	vare := vart*(q-vart);

	// Read in Crossingless matchings representation
	text := Read("Representation n="*IntegerToString(n)*" r="*IntegerToString(r));
	text := Split(text,"/");
	mats := [];

	// Get it to matrix
	for i:= 1 to #text do
	  mat := [];
	  for row in Split(text[i],"\n") do
	    out := Split(row,",");
	    newout := [];
	    for elt in out do
	      if elt eq "2" then
		newout := Append(newout,q+1);
	      else 
		if elt eq "1" then
		  newout := Append(newout,Sqrt(q));
		else
		  newout := Append(newout,q*0);
		end if;
	      end if;
	    end for;
	    mat := Append(mat,newout);
	  end for;
	  mats := Append(mats,mat);
	end for;

	// Correct the (1 + T_i) matrices to T_i matrices
	for j:=1 to #mats do
	  for i:=1 to #mats[j] do
	    mats[j][i][i] -:= 1;
	  end for;
	end for;


	// Read in Fibonacci representation and get it to matrix
	text := Read("Fib m="*IntegerToString(n+r)*" "*sp);
	text := Split(text,"/");
	fibmats := [];
	for i:= 1 to #text do
	  mat := [];
	  for row in Split(text[i],"\n") do
	    out := Split(row,",");
	    newout := [];
	    for elt in out do
	      case elt:
		when "a":
		  newout := Append(newout,vara); 
		when "b":
		  newout := Append(newout,varb); 
		when "c":
		  newout := Append(newout,varc); 
		when "d":
		  newout := Append(newout,vard); 
		when "e":
		  newout := Append(newout,vare); 
		when "0":
		  newout := Append(newout,L!0); 
	      end case;
	    end for;
	    mat := Append(mat,newout);
	  end for;
	  fibmats := Append(fibmats,mat);
	end for;

	// Convert to representation
	R := MatrixRing(L,#mats[1]);
	f := hom<H->R|mats>;
	V := RModule(L,#mats[1]);
	m := map<CartesianProduct(A,V)->V|t:->V!Vector(f(J(t[1]))*Transpose(Matrix(t[2\
	])))>; // the map AxV->V giving the module structure, basically just map to the matrix ring and then multiply with the vectors, plus coercing to different data types
	M:=Module(A,m);
	CompositionSeries(M);

	// print fibmats;

	// Get fibonacci representation
	Rfib := MatrixRing(L,#fibmats[1]);
	ffib := hom<H->Rfib|fibmats>;
	Vfib := RModule(L,#fibmats[1]);
	mfib := map<CartesianProduct(A,Vfib)->Vfib|t:->Vfib!Vector(ffib(J(t[1]))*Transpose(Matrix(t[2\
	])))>; // the map AxV->V giving the module structure, basically just map to the matrix ring and then multiply with the vectors, plus coercing to different data types
	Mfib := Module(A,mfib);

	return M,Mfib;
end function;

getIso := function(e,n,r,sp,H,A,T,q,L,J)
  M,MFib := getReps(e,n,r,sp,H,A,T,q,L,J);
  return IsIsomorphic(M,MFib);
end function;
      
e := 5;
n := 2;
s := 3;
r := s - n;
if r eq 0 or r eq 3 then
  sp := "star";
else
  sp := "p";
end if; 

H,A,T,q,L,J := generateHeckeAlgebra(e,s);
print "Algebra generated";

//M,MFib := getReps(e,n,r,sp,H,A,T,q,L,J);
print "n="*IntegerToString(n)*"\t r="*IntegerToString(r);
a,S := getIso(e,n,r,sp,H,A,T,q,L,J);
print S;
//if a then
//  print S;
//end if;

//for s := 4 to 7 do
//  H,A,T,q,L,J := generateHeckeAlgebra(e,s);
//  for n := s-3 to s do
//    if (n mod 2) eq 0 then
//      r := s - n;
//      if r eq 0 or r eq 3 then
//        sp := "star";
//      else
//        sp := "p";
//      end if; 
//      a,S := isIso(e,n,r,sp,H,A,T,q,L,J);
//      print "n="*IntegerToString(n)*"\t r="*IntegerToString(r);
//      if a then
//        print S;
//      end if;
//    end if;
//  end for;
//end for;
