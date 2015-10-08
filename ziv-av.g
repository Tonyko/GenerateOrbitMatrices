GraphFromAdjacencyMatrix := function (M)
  local n,i,j;
  n:=Size(M);
  return Graph(Group(()), [1..n], OnPoints, function(x,y) return M[x][y]<>0 ; end);
end;

GraphFromEdgeList := function(arg)
  local i,bl;
  if Size(arg)=2 and IsInt(arg[1]) then
    i:=arg[1];
    bl:=arg[2];
  else
    bl:=arg[1];
    i:=Maximum(Union(bl));
  fi ;
  return Graph(Group(()), [1..i], OnPoints, function(x,y) return [x,y] in bl; end);
end;


AsCoCo := function(a)
  if IsMatrix(a[1]) and DimensionsMat(a[1])[1]=DimensionsMat(a[1])[2] then
    return List(a, x -> GraphFromAdjacencyMatrix(x));
  fi;
  if IsList(a[1]) then
    return List(a, x -> GraphFromEdgeList(x));
  fi;
  if IsGraph(a[1]) then
    return a;
  fi;
end;


CMatrix := function(m,k)
  return List(m, x -> List(x, function(y) if y=k then return 1; fi; return 0; end));
end;


FromColorMatrix := function(m)
  local l,i;
  for i in [1..Size(m)] do
    m[i][i]:=m[i][i]-10000000000; 
  od ;
  l:=Union(m);
  return List(l, x -> CMatrix(m,x));
end;


CAut_Graph := function(c)
  local r,v,g,aut,t;
  
  aut:=Group(());
  if IsMatrix(c) then c:=FromColorMatrix(c); fi ;
  t:=AsCoCo(c);
  r:=List(t,DirectedEdges);
  v:=Concatenation([1..OrderGraph(t[1])],r,Tuples([1..OrderGraph(t[1])],2));
       
  g:=Graph( aut, v, function(x,y) 
    if IsInt(x) then return OnPoints(x,y) ; fi ;
    if IsInt(x[1]) then return OnTuples(x,y) ; fi ;
    return OnSetsTuples(x,y);
  end,
  function(x,y) 
    if IsInt(x) then return not IsInt(y) and x=y[1] ; fi ;
    if IsInt(y) then return not IsInt(x) and Size(x)=2 and y=x[2] ; fi ;
    return x in y or y in x ;
  end,true);
  return g;
end;


CAut_2 := function(c)
  return Action(AutomorphismGroup(CAut_Graph(c)),[1..Order(c)+Rank(c)]);
end;


IsIsomorphicCgr_2 := function(c1,c2)
  return IsIsomorphicGraph(CAut_Graph(c1),CAut_Graph(c2));
end;

