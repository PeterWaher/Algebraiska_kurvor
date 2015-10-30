
# FindGCD
FindGCD := proc(a, b)
	local a0 , b0 , m, n, r, GCD, BList, Position, A, B, switch;

	if not type(a, integer) then
		ERROR("a must be an integer!", a)
	end if;

	if not type(b, integer) then
		ERROR("b must be an integer!", b)
	end if;

	if a < b then
		m := b; 
		n := a; 
		switch := true
	else
		m := a; 
		n := b; 
		switch := false
	end if;

	a0 := m;
	b0 := n;
	GCD := n;
	r := m mod n;
	Position := 1;

	while r <> 0 do
		BList[Position] := -(m - r)/n;
		Position := Position + 1;
		m := n;
		n := r;
		GCD := r;
		r := m mod n
	end do;

	A := 0;
	B := 1;

	while 1 < Position do
		Position := Position - 1;
		r := B;
		B := A + BList[Position] * B;
		A := r
	end do;

	if switch then
		RETURN([GCD, B, A])
	else
		RETURN([GCD, A, B])
	end if
end proc;

# FindGCDList
FindGCDList := proc(IntegerList)
	local NrIntegers, Position, GCD, Coefficients, l;

	if not type(IntegerList, list) then
		ERROR("IntegerList must be a list of integers!", IntegerList)
	end if;

	NrIntegers := nops(IntegerList);
	if NrIntegers < 1 then
		ERROR("IntegerList must be a list of at least one integer!", IntegerList)
	elif NrIntegers = 1 then
		RETURN([IntegerList1, [1]])
	end if;

	l := FindGCD(IntegerList[1], IntegerList[2]);
	GCD := l[1];
	Coefficients := [l[2], l[3]];
	Position := 3;

	while Position <= NrIntegers do
		l := FindGCD(GCD, IntegerList[Position]);
		GCD := l[1];
		Coefficients := [op(Coefficients * l[2]), l[3]];
		Position := Position + 1
	end do;

	RETURN([GCD, Coefficients])
end proc;

# FindConductor
FindConductor := proc(Generators)
	local Conductor, l, MinValue, NrGenerators, g, ZMin, GCD, a, b, c, d, e, StartTime;

	StartTime := time();

	if not type(Generators, list) then
		ERROR("Generators must be a list of integers!", Generators)
	end if;

	l := FindGCDList(Generators);
	if l[1] <> 1 then
		ERROR("The Generators must not have a common divisor greater than 1!", [Generators, l[1]])
	end if;

	NrGenerators := nops(Generators);
	MinValue := min(op(Generators));
	ZMin := array(0..MinValue - 1);
	ZMin[0] := MinValue;

	for l from 1 to MinValue - 1 do
		ZMin[l] := 0
	end do;

	for l to NrGenerators do
		g := Generators[l];
		if g <> MinValue then
			GCD := gcd(g, MinValue);
			b := MinValue/GCD;
			for e from 0 to MinValue do
				if e = MinValue then
					c := 0
				elif ZMin[e] <> 0 then
					c := ZMin[e]
				else
					c := -1
				end if;

				if 0 <= c then
					for a to b do
						c := c + g;
						d := c mod MinValue;
						if ZMin[d] = 0 or c < ZMin[d] then
							ZMin[d] := c
						end if
					end do
				end if
			end do
		end if
	end do;

	Conductor := ZMin[0];
	for a to MinValue - 1 do
		if Conductor < ZMin[a] then
			Conductor := ZMin[a]
		end if
	end do;

	Conductor := Conductor - MinValue + 1;

	if nargs = 1 or args[2] then
		printf("Elapsed Time: %0.3f s.\n", time() - StartTime)
	end if;

	RETURN(Conductor)
end proc;

# FindSemiGroup
FindSemiGroup := proc(Generators)
	local Conductor, l, GCD, Generators2, MinValue, NrGenerators, g, ZMin, GCD2,a, b, c, d, e, SemiGroup, StartTime;
	StartTime := time();

	if not type(Generators, list) then
		ERROR("Generators must be a list of integers!", Generators)
	end if;

	l := FindGCDList(Generators);
	GCD := l[1];
	if GCD = 1 then
		Generators2 := Generators
	else
		Generators2 := Generators/GCD
	end if;

	NrGenerators := nops(Generators2 );
	MinValue := min(op(Generators2 ));
	ZMin := array(0..MinValue - 1);
	ZMin[0] := MinValue;

	for l from 1 to MinValue - 1 do
		ZMin[l] := 0
	end do;

	for l to NrGenerators do
		g := Generators2[l];
		if g <> MinValue then
			GCD2 := gcd(g, MinValue);
			b := MinValue/GCD2;

			for e from 0 to MinValue do
				if e = MinValue then
					c := 0
				elif ZMin[e] <> 0 then
					c := ZMin[e]
				else
					c := -1
				end if;

				if 0 <= c then
					for a to b do
						c := c + g;
						d := c mod MinValue;

						if ZMin[d] = 0 or c < ZMin[d] then
							ZMin[d] := c
						end if
					end do
				end if
			end do
		end if
	end do;

	Conductor := ZMin[0];
	for a to MinValue - 1 do
		if Conductor < ZMin[a] then
			Conductor := ZMin[a]
		end if
	end do;

	Conductor := Conductor - MinValue + 1;
	ZMin := array(1..Conductor);
	for l to Conductor do
		ZMin[l] := 0
	end do;

	for l to NrGenerators do
		g := Generators2[l];

		for a to Conductor + 1 do
			if Conductor < a then
				b := g
			elif ZMin[a] <> 0 then
				b := ZMin[a] + g
			else
				b := 0
			end if;

			if 0 < b then
				while b <= Conductor do
					ZMin[b] := b;
					b := b + g
				end do
			end if
		end do
	end do;

	SemiGroup := {};

	for a to Conductor do
		if 0 < ZMin[a] then
			SemiGroup := SemiGroup union {ZMin[a]}
		end if
	end do;

	if GCD <> 1 then
		SemiGroup := GCD * SemiGroup
	end if;

	if nargs = 1 or args2 then
		printf("Elapsed Time: %0.3f s.\n", time() - StartTime)
	end if;

	RETURN([Conductor * GCD, SemiGroup])
end proc;

# FindSemiGroupFromPolynomialRing
FindSemiGroupFromPolynomialRing := proc(PolynomialGenerators, Variable)
	local StartTime, GCD, Conductor, Polynomials, NrPolynomials, 
		Generators, PolynomialsByOrder, Orders, ExplicitNotation,
		FirstExplicitNotation, NrOriginalPolynomials, OriginalOrders,
		MaxExponent, MinOrder, HighestOrder, CalcExplicit, MaxTerm,
		MaxTerms, a, b, c, d, e, f, g, h, e1, e2, g1, g2, g3, d1, d2, d3, p;

	StartTime := time();

	CalcExplicit := 4 <= nargs and args[4];

	if not type(PolynomialGenerators, list) then
		ERROR("PolynomialGenerators must be a list of polynomials!", PolynomialGenerators)
	end if;

	NrOriginalPolynomials := nops(PolynomialGenerators);
	if NrOriginalPolynomials = 0 then
		ERROR("List of polynomials cannot be empty!")
	end if;

	f := PolynomialGenerators;
	for a to NrOriginalPolynomials do
		g1 := f[a];
		if not type(g1, polynom) then
			ERROR("PolynomialGenerators must be a list of polynomials!", PolynomialGenerators)
		end if;

		ExplicitNotation[a] := p[a];
		OriginalOrders[a] := ldegree(g1 , Variable);
		Orders[a] := OriginalOrders[a]
	end do;

	for a from 2 to NrOriginalPolynomials do
		for b to a - 1 do
			if Orders[a] < Orders[b] then
				c := Orders[a];
				Orders[a] := Orders[b];
				Orders[b] := c;
				c := f[a];
				f[a] := f[b];
				f[b] := c;
				c := ExplicitNotation[a];
				ExplicitNotation[a] := ExplicitNotation[b];
				ExplicitNotation[b] := c
			end if
		end do
	end do;

	MinOrder := Orders[1];
	GCD := MinOrder;
	Generators := [MinOrder];
	g1 := f[1];
	e := coeff(g1, Variable, MinOrder);
	Polynomials[1] := g1/e;
	ExplicitNotation[1] := ExplicitNotation[1]/e;
	g := 1;

	for a from 2 to NrOriginalPolynomials do
		g1 := f[a];
		b := ExplicitNotation[a];

		for h to g do
			g2 := Polynomials[h];
			d2 := ldegree(g2 , Variable);
			e := coeff(g1, Variable, d2);

			if e <> 0 then
				g1 := g1 - e * g2;
				b := b - e * ExplicitNotation[h]
			end if
		end do;

		if g1 <> 0 then
			d := ldegree(g1 , Variable);
			Generators := [op(Generators), d];
			GCD := gcd(GCD, d);
			e := coeff(g1 , Variable, d);
			g1 := g1/e;
			b := b/e;
			g := g + 1;
			Polynomials[g] := g1;
			ExplicitNotation[g] := b;

			for h to g - 1 do
				g2 := Polynomials[h];
				e := coeff(g2, Variable, d);
				if e <> 0 then
					g2 := g2 - e * g1;
					Polynomials[h] := g2;
					ExplicitNotation[h] := ExplicitNotation[h] - e * b
				end if
			end do
		end if
	end do;

	NrPolynomials := g;
	if GCD = 1 then
		Conductor := FindConductor(Generators, false);
		for c to NrOriginalPolynomials do
			MaxExponent[c] := ceil((Conductor + 1)/OriginalOrders[c])
		end do
	else
		Conductor := infinity
	end if;

	HighestOrder := max(op(Generators));
	for a from MinOrder to HighestOrder do
		PolynomialsByOrder[a] := 0
	end do;

	for a to NrPolynomials do
		f := Polynomials[a];
		d := ldegree(f, Variable);
		Orders[a] := d;
		PolynomialsByOrder[d] := a;
		FirstExplicitNotation[d] := ExplicitNotation[a]
	end do;

	a := 1;
	while a <= NrPolynomials do
		g1 := Polynomials[a];
		d1 := Orders[a];

		if g1 <> 0 and d1 <= Conductor then
			if CalcExplicit then
				e1 := ExplicitNotation[a]
			end if;

			b := 1;
			while b <= a do
				if b = a then
					g2 := g1;
					d2 := d1
				else
					g2 := Polynomials[b];
					d2 := Orders[b]
				end if;

				if g2 <> 0 and d2 <= Conductor then
					d := sort(simplify(expand(g1 * g2)));
					if CalcExplicit then
						e2 := expand(e1 * ExplicitNotation[b])
					end if;

					d3 := d1 + d2;
					while d3 <= HighestOrder and d3 <= Conductor and d <> 0 and PolynomialsByOrder[d3] <> 0 do
						c := PolynomialsByOrder[d3];
						g3 := Polynomials[c];

						if g3 = 0 then
							ERROR("Runtime Error", d)
						end if;

						e := coeff(d, Variable, d3);
						d := d - e * g3;

						if CalcExplicit then
							e2 := e2 - e * ExplicitNotation[c]
						end if;

						d3 := ldegree(d, Variable)
					end do;

					if d <> 0 and d3 <= Conductor then
						e := coeff(d, Variable, d3);
						d := d/e;

						if CalcExplicit then
							e2 := e2/e
						end if;

						Generators := [op(Generators), d3];
						GCD := gcd(GCD, d3);

						if GCD = 1 and Conductor = infinity then
							Conductor := FindConductor(Generators, false);

							if Conductor < d3 then
								Conductor := d3
							end if;

							for c to NrOriginalPolynomials do
								MaxExponent[c] := ceil((Conductor + 1)/OriginalOrders[c]);
								MaxTerms[c] := p[c]^MaxExponent[c]
							end do;

							MaxTerm := Variable^(Conductor+1);
							if Conductor < degree(d, Variable) then
								d := rem(d, MaxTerm, Variable);
								if CalcExplicit then
									for c to NrOriginalPolynomials do
										if MaxExponent[c] < degree(e2 , p[c]) then
											e2 := expand(rem(e2 , MaxTerms[c], p[c]))
										end if
									end do
								end if
							end if;

							for c to NrPolynomials do
								e := Polynomials[c];
								if e <> 0 and Conductor < degree(e, Variable) then
									Polynomials[c] := rem(e, MaxTerm, Variable);

									if CalcExplicit then
										if Polynomials[c] = 0 then
											ExplicitNotation[c] := 0
										end if
									else
										e := ExplicitNotation[c];
										for f to NrOriginalPolynomials do
											if MaxExponent[f] < degree(e, p[f]) then
												e := expand(rem(e, MaxTerms[f] , p[f]))
											end if
										end do;

										ExplicitNotation[c] := e
									end if
								end if
							end do
						elif Conductor < degree(d, Variable) then
							d := rem(d, MaxTerm, Variable);
							if CalcExplicit then
								for c to NrOriginalPolynomials do
									if MaxExponent[c] < degree(e2 , p[c]) then
										e2 := expand(rem(e2 , MaxTerms[c], p[c]))
									end if
								end do
							end if
						end if;

						NrPolynomials := NrPolynomials + 1;
						Polynomials[NrPolynomials] := d;
						Orders[NrPolynomials] := d3;

						if CalcExplicit then
							ExplicitNotation[NrPolynomials] := e2;
							FirstExplicitNotation[d3] := e2
						end if;

						if HighestOrder < d3 then
							HighestOrder := HighestOrder + 1;
							while HighestOrder < d3 do
								PolynomialsByOrder[HighestOrder] := 0;
								HighestOrder := HighestOrder + 1
							end do
						end if;

						PolynomialsByOrder[d3] := NrPolynomials;
						h := NrPolynomials - 1;

						for c to h do
							e := Polynomials[c];
							if e <> 0 then
								f := coeff(e, Variable, d3);
								if f <> 0 then
									e := e - f * d;
									NrPolynomials := NrPolynomials + 1;
									g := Orders[c];
									Polynomials[NrPolynomials] := e;
									Orders[NrPolynomials] := g;
									PolynomialsByOrder[g] := NrPolynomials;

									if CalcExplicit then
										e := ExplicitNotation[c] - f * e2;
										ExplicitNotation[NrPolynomials] := e
									end if;

									Polynomials[c] := 0
								end if
							end if
						end do
					end if
				end if;

				b := b + 1
			end do
		end if;

		a := a + 1
	end do;

	if Conductor = infinity then
		ERROR("No upper bound found.", PolynomialGenerators)
	end if;

	Generators := sort(Generators);
	a := {};
	b := [];
	for c to nops(Generators) do
		d := Generators[c];
		if d <= Conductor and not member(d, a) then
			b := [op(b), d];
			e := d;
			f := {};
			while e <= Conductor do
				f := f union {e};

				for h to nops(a) do
					g := a[h] + e;
					if g <= Conductor then 
						f := f union {g}
					end if
				end do;

				e := e + d
			end do;

			a := a union f
		end if
	end do;

	Generators := b;
	if CalcExplicit then
		for c to NrOriginalPolynomials do
			MaxExponent[c] := ceil((Conductor + 1)/OriginalOrders[c]);
			MaxTerms[c] := p[c]^MaxExponent[c]
		end do;

		MaxTerm := Variable^(Conductor +1);

		for a to nops(Generators) do
			b := Generators[a];
			c := FirstExplicitNotation[b];

			for e to NrOriginalPolynomials do
				c := expand(rem(c, MaxTerms[e], p[e]))
			end do;

			d := convert(c, list);
			if 1 < nops(d) then
				for e to nops(d) do
					f := d[e];

					for g to NrOriginalPolynomials do
						f := subs(p[g] = PolynomialGenerators[g], f)
					end do;

					if b < ldegree(f) then 
						c := c - d[e] 
					end if
				end do
			end if;

			d := c;
			for g to NrOriginalPolynomials do
				d := subs(p[g] = PolynomialGenerators[g], d)
			end do;

			d := sort(expand(d));
			e := lcm(denom(c), denom(d));
			c := sort(e * c);
			d := sort(e * d);
			print(c = d)
		end do
	end if;

	if nargs < 3 or args[3] then
		printf("Elapsed Time: %0.3f s.\n", time() - StartTime)
	end if;

	RETURN(Generators)
end proc;
