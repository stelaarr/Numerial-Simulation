### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 7e95bb80-a714-11ed-29d5-0bcfe8e0bbbf
md"""###
Ejercicio 4

(a) Implementa el algoritmo general para el cálculo de los coeficientes de una fórmula de diferencias finitas para aproximar la derivada de orden $k$ de $u$ en el punto $x$, para un conjunto de nodos xnod de tamaño $n>k$. Para ello, define una función: $fdcoeff(k,x,xnod)$ que devuelva los coeficientes de la fórmula de diferencias finitas.

(b) Comprueba que la solución es idéntica si se evalúa la derivada en cero y se trasladan los nodos al origen, es decir, comprueba que:

$fdcoeff(k,0,xnod-x) = fdcoeff(k,x,xnod)$

de modo que podemos prescindir de la variable $x$."""

# ╔═╡ 98393f40-654e-4786-960f-89383b0fc610
md"""### Resolucion:"""

# ╔═╡ f2dd77fd-4556-4a6d-af34-b6bf9c198a58
function fdcoeff(k,x,xnod)
	n=length(xnod);
	A=zeros(n,n);
	B=zeros(n,1);
	B[k+1,1]=1
	for i in 1:n
		A[1,i]=1.;
	end
	for row in 2:n
		for col in 1:n
			A[row,col]=(1/factorial(row-1))*((xnod[col]-x)^(row-1))
		end
	end
	coeffi=inv(A)*B
	return coeffi
end
		

# ╔═╡ 356ef9fd-0295-40b9-b5d5-9bc610966314
# una prueba para calcular los coeficientes de la derivada segunda centrada
fdcoeff(2,0,[-1,0.0,1])

# ╔═╡ a72bf2e0-8754-445c-a6a2-6dfc5662ad47
# la derivada primera con el esquema de tercer orden
fdcoeff(1,0,[-2,-1,0,1])

# ╔═╡ 854c4d34-c35d-49c1-b830-34317f6e6aa1
#prueba para ver si desplazada da el mismo resultado
fdcoeff(1,1,[-1,0,1,2])

# ╔═╡ Cell order:
# ╠═7e95bb80-a714-11ed-29d5-0bcfe8e0bbbf
# ╠═98393f40-654e-4786-960f-89383b0fc610
# ╠═f2dd77fd-4556-4a6d-af34-b6bf9c198a58
# ╠═356ef9fd-0295-40b9-b5d5-9bc610966314
# ╠═a72bf2e0-8754-445c-a6a2-6dfc5662ad47
# ╠═854c4d34-c35d-49c1-b830-34317f6e6aa1
