
# MakePoly

MakePoly := proc(a, Variable)
	local i, r;

	r := 0;
	for i to nops(a) do 
		if op(i, a) <> 0 then
			r := r + op(i, a) * Variable^(i-1)
		fi
	od;

	RETURN(sort(r))
end;

# Composition

Composition := proc(f, g, Variable, MaxDegree)
	local gc, gcn, gcnt, i, j, k, r, c, floop;

	gc := array(0..MaxDegree);
	gcn := array(0..MaxDegree);
	gcnt := array(0..MaxDegree);
	r := array(0..MaxDegree);

	for i from 0 to MaxDegree do
		gc[i] := coeff(g, Variable, i); 
		gcn[i] := gc[i]; 
		r[i] := 0
	od;

	floop := degree(f);
	if MaxDegree < floop then 
		floop := MaxDegree
	fi;

	r[0] := coeff(f, Variable, 0);
	for i to floop do
		c := coeff(f, Variable, i);
		for j from 0 to MaxDegree do 
			r[j] := r[j] + c * gcn[j] 
		od;

		if i < MaxDegree then
			for j from 0 to MaxDegree do
				gcnt[j] := 0; 
				for k from 0 to j do
					gcnt[j] := gcnt[j] + gc[k] * gcn[j-k]
				od
			od;

			for j from 0 to MaxDegree do 
				gcn[j] := gcnt[j] 
			od
		fi
	od;

	RETURN(eval(r))
end;

# Reparametrize

Reparametrize := proc(x, y, Variable, t0 , MaxDegree, Branch, AllowNegation)
	local xt, yt, a, p, q, x0 , y0 , i, j, k, shifted, temp, sols, StartTime, MinDegree,
		arg, z, Negate;

	StartTime := time();
	if not type(Branch, integer) then
		ERROR("Branch must be an integer!", Branch)
	fi;

	if t0 <> 0 then
		xt := collect(expand(convert(taylor(
			collect(expand(subs(Variable = Variable + t0, x)), Variable),
			Variable = 0, MaxDegree + 1), polynom)), Variable);

		yt := collect(expand(convert(taylor(
			collect(expand(subs(Variable = Variable + t0, y)), Variable),
			Variable = 0, MaxDegree + 1), polynom)), Variable)
	else
		xt := collect(expand(convert(taylor(x, Variable = 0, MaxDegree + 1),
			polynom)), Variable);

		yt := collect(expand(convert(taylor(y, Variable = 0, MaxDegree + 1),
			polynom)), Variable);
	fi;

	x0 := coeff(xt, Variable, 0);
	y0 := coeff(yt, Variable, 0);
	xt := xt - x0;
	yt := yt - y0;
	MinDegree := ldegree(xt, Variable);
	shifted := ldegree(yt, Variable) < MinDegree;

	if shifted then
		temp := xt; 
		xt := yt; 
		yt := temp; 
		MinDegree := ldegree(xt, Variable)
	fi;

	a := array(1..MaxDegree - MinDegree + 1);
	z := 1/coeff(xt, Variable, MinDegree);
	if z < 0 and AllowNegation then 
		xt := -xt; 
		z := -z; 
		Negate := true
	else 
		Negate := false
	fi;

	arg := argument(z);
	a[1] := abs(z)^(1/MinDegree) * (
		cos((2 * (MinDegree - 2 + Branch) * Pi + arg)/MinDegree)
		+ I * sin((2 * (MinDegree - 2 + Branch) * Pi + arg)/MinDegree));

	p := a[1] * Variable;
	for j from 2 to MaxDegree - MinDegree + 1 do 
		p := p + a[j] * Variable^j 
	od;

	q := Composition(yt, p, Variable, MaxDegree);
	p := Composition(xt, p, Variable, MaxDegree);
	i := MinDegree + 1;
	j := 2;

	while i <= MaxDegree do
		sols := [solve(p[i] = 0, a[j])];
		q[i] := expand(subs(a[j] = sols[1], q[i]));

		if i < MaxDegree then
			for k from i + 1 to MaxDegree do
				p[k] := subs(a[j] = sols[1], p[k]); 
				q[k] := subs(a[j] = sols[1], q[k])
			od
		fi;

		a[j] := sols[1];
		i := i + 1;
		j := j + 1
	od;

	xt := q[0];
	for i to MaxDegree do 
		if q[i] <> 0 then 
			xt := xt + q[i] * Variable^i 
		fi 
	od;

	if shifted then
		xt := x0 + xt;
		if Negate then 
			yt := y0 - Variable^MinDegree
		else 
			yt := y0 + Variable^MinDegree
		fi
	else
		yt := y0 + xt;
		if Negate then 
			xt := x0 - Variable^MinDegree
		else 
			xt := x0 + Variable^MinDegree
		fi
	fi;

	printf("Elapsed Time: %0.3f s.\n", time() - StartTime);

	RETURN([xt, yt])
end;
