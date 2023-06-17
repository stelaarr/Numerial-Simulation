### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 159bb200-969e-4928-a6d8-9e43fe7d07ac
using Plots

# ╔═╡ cf2dfedc-f412-4261-a5d9-9e557eed5c42
using LinearAlgebra

# ╔═╡ 6892b420-b0fe-11ed-1458-711f19c46668
md"""
# Practica 3: Problemas estacionarios II
"""

# ╔═╡ eb9c4ea1-308c-424c-8b1d-7b936f68b66e
md""" ALGORITMO DE THOMAS PARA MATRICES TRIDIAGONALES"""

# ╔═╡ a39b572a-3e38-4b10-b215-1855555b1cd4
function tridiag(a,b,c,y)
	n = length(y)
	x = zeros(n)
	bg = zeros(n)
	yg = zeros(n)	
	bg[1] = b[1]
	yg[1] = y[1]
	for k = 1:n-1
		bg[k+1] = b[k+1] - (a[k+1]*c[k])/bg[k]
		yg[k+1] = y[k+1] - (a[k+1]*yg[k])/bg[k]
	end
	x[n] = yg[n]/bg[n]
	for k = n-1:-1:1
		x[k] = (yg[k]-c[k]*x[k+1])/bg[k]
	end
	return x
end

# ╔═╡ d4de887a-615a-4d38-8b62-5946b27a9185
md""" # Ejercicio 1"""

# ╔═╡ 743a19df-7b68-48bb-abe2-1480a445e7e8
md""" Considera el siguiente problema con condiciones de contorno mixtas:

$u''(x)=e^{x},\;\;u'(0)=0,\;\;u(1)=3$

La solución exacta viene dada por:

$u(x) = e^{x} - e + 4 - x$

Para resolver el problema, utiliza el esquema general de segundo orden para las derivadas y los siguientes esquemas para la condición de Neumann. Estima el orden del método en la norma del máximo en cada caso con distintos valores de $h$ como en la hoja anterior.

(a) Esquema de primer orden para la derivada en $x=0$:

$\frac{U_{1} - U_{0}}{h} = \alpha$

El sistema resultante tiene $m+2$ ecuaciones, que pueden reducirse usando las condiciones de contorno."""

# ╔═╡ 266b104c-c0ec-4fa5-b706-7d2566f3c531
md""" # Solución:"""

# ╔═╡ 28dd18e2-6651-4e32-9af0-682c5ce9d521
md""" Solución exacta"""

# ╔═╡ d5345ec5-2771-49a6-83dc-79e598493460
function solexact(x)
	n = exp(x) - exp(1) + 4 - x
	return n
end

# ╔═╡ e0050aa5-5403-418f-b44e-e9bd88b51f88
md""" Función para Condiciones de Neumann izquierda con método 1:"""

# ╔═╡ b5520aa5-63ec-41ff-9f4f-20ddd85d5e9c
function ECCN_met1(f, alpha, beta, m, x)
	#defino primero el paso y los nodos
	h=1/(m+1);
	x=zeros(m+2);
	x=0:h:1;
	#defino ahora las diagonales
	a=ones(m)/h^2;
	b=(-2)*ones(m)/h^2;
	b[1]=-1/h^2; 
	c=ones(m)/h^2;
	#defino el vector y
	y=zeros(m);
	y=f.(x[2:1:end-1]); # el f. aplica a cada componente del array
	y[1]+= alpha/h;
	y[end]-= beta/h^2;
	#Entonces la solucion u será
	u=tridiag(a,b,c,y)
	#quiero que me devuleva los nodos x y la solución
	return x, [u[1] - alpha*h;u;beta]
end

# ╔═╡ a206d82e-010c-486a-87be-f90bd4c114ef
md""" Tengo que probar con $h=1/10,1/20,1/40$"""

# ╔═╡ 10d5d1e2-e085-4acb-8e90-50799cb91fe0
xx,soll=ECCN_met1(exp,0,3,9,0)

# ╔═╡ 4a7679e7-edf6-4ad3-8397-326c05f683e4
xx1,soll1=ECCN_met1(exp,0,3,19,0);

# ╔═╡ 1cd9579e-59b8-49ac-886b-3b24b1faf030
xx2,soll2=ECCN_met1(exp,0,3,39,0);

# ╔═╡ 49ff9dbe-e0da-4072-8dcd-c693feaad5f0
begin
	nn=zeros(length(xx));
	nn=solexact.(xx);
end

# ╔═╡ 42724883-c2b4-4d3b-b2c3-b9819be1b5aa
plot(xx,[soll,nn],title = "Con 10 nodos internos", label =["Solución numérica U" "Solución exacta u"])

# ╔═╡ f61fe6d6-7724-4368-bcc7-1c28269318f3
md""" Recordando como calcular el orden en las notas"""

# ╔═╡ 5690b3cb-8fbc-46f3-9a4e-e8e7be1b8919
begin
	uh=soll
	uh2=soll1[1:2:length(soll1)]
	uh4=soll2[1:4:length(soll2)]
	norma2_uh=norm(uh-uh2)
	norma2_uh2=norm(uh2-uh4)
end

# ╔═╡ f29616f9-459f-4ddc-803c-04eea8d41f48
log2(norma2_uh/norma2_uh2)

# ╔═╡ fad3e8ba-5797-42a4-aa13-59d7f5033eb0
md"""(b) Esquema de segundo orden, aplicando la derivada segunda en 
 y la condición de contorno.

El sistema resultante tiene $m+3$ ecuaciones que pueden reducirse eliminando $U_{-1}$ entre las dos ecuaciones anteriores."""

# ╔═╡ ac812159-7bc1-4fea-bd1d-323edc85b1cf
md""" # Solución:"""

# ╔═╡ ba23629d-0664-4af3-aba0-cebcc2e332a9
md""" Defino función para el metodo 2:"""

# ╔═╡ 41246027-3220-4256-bd2e-f8c5f34c8176
function ECCN_met2(f,alpha,beta,m,x)
	h = 1/(m + 1)
	x = zeros(m+2)
	x = 0:h:1
	
    a = ones(m+1)/h^2
    b = (-2)*ones(m+1)/h^2
	b[1] = -1/h^2
	c = ones(m+1)/h^2
    y = zeros(m+1)

	y = f.(x[2:end-1])
	y[1] += 1/2 + alpha/h
    y[end] -= beta/h^2

	u = tridiag(a,b,c,y)

    return x,[u[1]- alpha*h - y[1]*h^2/2;u;beta]
end

# ╔═╡ bdc568c3-16ee-4b9f-ad42-a53067e8ebe2
xm2,solm2=ECCN_met2(exp,0,3,9,0)

# ╔═╡ b2b6d7dc-b9d1-4b94-8af0-624d59e72e17
xxm2,sollm2=ECCN_met2(exp,0,3,19,0)

# ╔═╡ 704af3a8-5836-4400-8cd9-6b30efa514ee
xxxm2,solllm2=ECCN_met2(exp,0,3,39,0)

# ╔═╡ 35faa2fd-e37c-4f9a-87cb-90fe2513d049
plot(xm2,[solm2,nn], title = "Con 10 nodos internos", label =["Solución numérica U" "Solución exacta u"])

# ╔═╡ 9b01627c-b6a4-4865-96b0-9f6483616eb4
md""" Ahora el orden será: (Defino función orden que es más útil)"""

# ╔═╡ 5b7e8797-2b96-4ee7-ae72-4a333369d289
function order(sol1,sol2,sol3)
	e1 = sol1 - sol2[1:2:length(sol2)];
	e2 = sol2[1:2:length(sol2)] -  sol3[1:4:length(sol3)];
	norma1 = norm(e1);
	norma2 = norm(e2);
	p = log2(norma1/norma2)
	return p
end

# ╔═╡ 8d9dd22f-d8cd-4720-8248-e4b58dd47a83
order(solm2,sollm2,solllm2)

# ╔═╡ acd29996-cfc4-4e16-baae-5ffc6ea26a47
md""" (c) Esquema a un lado para la derivada primera en $x=0$:"""

# ╔═╡ 7bd417ef-db85-472d-bcd9-0f2a5f8f43ba
md""" # Solución:"""

# ╔═╡ 70f3183b-2818-4249-8414-6a396f0d073f
md""" Defino función para el tercer método"""

# ╔═╡ ed4365e0-205a-4f9d-b2a0-cb22ba65c730
function ECCN_met3(f, alpha, beta, m, x)
	#defino primero el paso y los nodos
	h=1/(m+1);
	x=zeros(m+2);
	x=0:h:1;
	#defino ahora las diagonales
	a=(1)*ones(m)/h^2;
	b=(-2)*ones(m)/h^2;
	b[1]=-2/(3*(h^2)); 
	c=(1)*ones(m)/h^2;
	c[1]=2/(3*(h^2)); 
	#defino el vector y
	y=zeros(m+1);
	y=f.(x[2:end-1]); # el f. aplica a cada componente del array
	y[1]+= (2*alpha)/(3*h);
	y[end]-= beta/h^2;
	#Entonces la solucion u será
	u=tridiag(a,b,c,y)
	#quiero que me devuleva los nodos x y la solución
	return x, [(4*u[1] - 2*alpha*h - u[2])/3;u;beta]
end

# ╔═╡ 3503d295-7f6b-4f9b-848b-6159317fd19a
xm3,solm3=ECCN_met3(exp,0,3,9,0);

# ╔═╡ c1741f3d-9e66-4b10-9e37-d3606842f704
xxm3,sollm3=ECCN_met3(exp,0,3,19,0);

# ╔═╡ 6b87c10e-1856-47ff-994a-d83bc6d01ace
xxxm3,solllm3=ECCN_met3(exp,0,3,39,0);

# ╔═╡ 67c329be-ea7f-4dde-bd80-7fe556f7a74a
plot(xm3,[soll,solm2,solm3,nn], title = "Con 10 nodos internos", label =["Solución (a)" "Solución (b)" "Solución (c)" "Solución exacta u"])

# ╔═╡ fd007c48-34fc-4fee-b936-f66ce9f6128e
md""" El orden es:"""

# ╔═╡ f9ad13d1-0645-4482-813a-dfdcd16ac247
order(solm3,sollm3,solllm3)

# ╔═╡ 43b1bb5f-a431-4a7f-980c-150ea0db6473
md""" # Ejercicio 2"""

# ╔═╡ cfa02242-2d07-4e40-aea5-d2dd8f5032dd
md""" Implementa el esquema general para una ecuación lineal de segundo orden para el caso de la ecuación de la propagación del calor con condiciones de Dirichlet:

$(\kappa(x)u'(x))' = -sin(\pi x),\;\;u(0)=2,\;\;u(1)=1$

donde 

$\kappa(x)=1-(x-\frac{1}{2})^2$

(a) Implementa el esquema centrado, resolviendo el sistema lineal para distintos valores de $h$.Estima el orden del método."""

# ╔═╡ ebc4c77e-025c-416a-9b0d-f4e7de796928
md""" # Solución:"""

# ╔═╡ 9de8ba60-db2f-41a7-a31d-7579f5158bd0
md""" Resolviendo queda:

$\kappa(x)u''(x) + k'(x)u'(x) = -sin(\pi x)$"""

# ╔═╡ 23e684b0-f03d-4a92-aee2-4d63960ba796
md""" Defino función para resolver ecuación general con Condiciones de Dirichlet:"""

# ╔═╡ 8d938768-d418-4b89-a9b7-d27c9f748730
function EGCD(lambda,mu,eta,f,intervalo,alpha,beta,m)
	#defino primero el paso y los nodos, las diagonales y el vector  y
	h = (intervalo[2]-intervalo[1])/(m+1)
    a = ones(m)
    b = ones(m)
    c = ones(m)
    y = zeros(m)
	x = zeros(m)
	for i = 1:m
		x[i] = intervalo[1]  + i*h
		a[i] = lambda(x[i])/h^2 - mu(x[i])/(2*h)
		b[i] = (-2)*lambda(x[i])/h^2 + eta(x[i])
		c[i] = lambda(x[i])/h^2 + mu(x[i])/(2*h)
		y[i] = f.(x[i])
	end
	y[1] -= (lambda(x[1])/h^2 - mu(x[1])/(2*h))*alpha
	y[end] -= (lambda(x[end])/h^2 + mu(x[end])/(2*h))*beta
	u = tridiag(a,b,c,y)
	return [intervalo[1];x;intervalo[2]],[alpha;u;beta]
end

# ╔═╡ 8dc14c5e-b6a3-480a-9d30-3476805557ee
xex2,sol_ex2=EGCD(x->1+(x-1/2)^2,x->2*(x-1/2),x->0,x->-sin(pi*x),[0,1],2,1,9);

# ╔═╡ b69d1fe6-445c-4fb8-9a4c-a245799b3e28
xxex2,soll_ex2=EGCD(x->1+(x-1/2)^2,x->2*(x-1/2),x->0,x->-sin(pi*x),[0,1],2,1,19);

# ╔═╡ cb05c42b-afb0-4c73-aac1-8ebfcacc064b
xxxex2,solll_ex2=EGCD(x->1+(x-1/2)^2,x->2*(x-1/2),x->0,x->-sin(pi*x),[0,1],2,1,39);

# ╔═╡ e1810634-1ff5-43c4-bf81-5e6870c07a25
md""" El orden será:"""

# ╔═╡ 3da56891-d0ce-4287-9e80-4105a6da99ff
order(sol_ex2,soll_ex2,solll_ex2)

# ╔═╡ 379238de-646c-4db9-a853-7593616bd1c9
md""" (b) Implementa el esquema conservativo basado en volúmenes finitos. Estima el orden del método."""

# ╔═╡ 4ad30e26-755d-4cf1-a5f3-af714580ab5f
md""" # Solución:"""

# ╔═╡ fe5cb681-30c8-4dca-8744-955a80cb1058
md""" Metodo de volumenes finitos"""

# ╔═╡ ffcb9628-d626-42c2-8519-d9df668663df
function EGVF(k,f,intervalo,alpha,beta,m)
	#defino primero el paso y los nodos, las diagonales y el vector  y
	h = (intervalo[2]-intervalo[1])/(m+1)
    a = ones(m)
    b = ones(m)
    c = ones(m)
    y = zeros(m)
	x = (intervalo[1] + h/2):h:(intervalo[2] + h/2)
	
	a = k.(x[1:end-2])/h^2 
	b = - (k.(x[1:end-2]) + k.(x[2:end-1]))/h^2
	c = k.(x[2:end-1])/h^2
	y = (f.(x[1:end-2])+f.(x[2:end-1]))/2
	y[1] -= (k(x[1])/h^2)*alpha
	y[end] -= (k(x[end-1])/h^2)*beta
	u = tridiag(a,b,c,y)
	return x,[alpha;u;beta]
end

# ╔═╡ e9068806-df9f-4990-987a-792a3fa9d34d
xc,solc=EGVF(x->1+(x-1/2)^2,x->-sin(pi*x),[0,1],2,1,9);

# ╔═╡ f2f75266-6124-4656-bdbe-8bdaecb02cf0
xxc,sollc=EGVF(x->1+(x-1/2)^2,x->-sin(pi*x),[0,1],2,1,19);

# ╔═╡ 5240a2a5-2da0-4c66-9d48-c9baf3762e29
xxxc,solllc=EGVF(x->1+(x-1/2)^2,x->-sin(pi*x),[0,1],2,1,39);

# ╔═╡ b8376f9f-a95a-4bb8-bb00-fd149649be91
plot(xc,[sol_ex2 solc], label = ["Esquema centrado 2do orden" "Esquema conservativo"])

# ╔═╡ 03539f0a-b958-49d9-a80e-8c94ebae688e
md""" El orden será:"""

# ╔═╡ da37f334-8b65-46cd-b085-9dca5bfe6d3f
order(solc,sollc,solllc)

# ╔═╡ 7859444f-9afb-4e24-8909-e6a2665e5242
md"""(c) Resuelve el problema con condición de Neumann a la izquierda: $u'(0)=0$"""

# ╔═╡ e3e5ff04-8683-4410-8204-bda2807a9e8b
md""" # Solución:"""

# ╔═╡ 46bbccc1-1a57-419c-b30d-30e24bea8cf4
function EGVF_CD(k,f,intervalo,alpha,beta,m)
	#defino primero el paso y los nodos, las diagonales y el vector  y
	h = (intervalo[2]-intervalo[1])/m
    a = ones(m)
    b = ones(m)
    c = ones(m)
    y = zeros(m)
	x = (intervalo[1] + h/2):h:(intervalo[2] + h/2)
	
	a = k.(x[1:end-2])/h^2 
	b = - (k.(x[1:end-2]) + k.(x[2:end-1]))/h^2
	b[1]= - k(x[2])/h^2
	b[end]= - (2*k(x[end-1]) + k(x[end-2]))/h^2
	c = k.(x[2:end-1])/h^2
	c[1]=  k(x[2])/h^2
		
	y = (f.(x[1:end-2])+f.(x[2:end-1]))/2
	y[1] += (k(x[1])/h^2)*alpha
	y[end] -= (2*k(x[end-1])/h^2)*beta
	u = tridiag(a,b,c,y)
	return x[1:end-2],u
end

# ╔═╡ 6bed92ed-1234-4081-b386-29cc7eb1281a
x1,y1=EGVF_CD(x->1+(x-1/2)^2,x->-sin(pi*x),[0,1],0,1,9)

# ╔═╡ 3e4c1677-bee6-485e-bd20-08d68bb7972f
plot(x1,y1)

# ╔═╡ 02f1147a-7995-429b-a56d-564548e649f9
md""" # Ejercicio 3:"""

# ╔═╡ 7319f969-ae44-4363-8b65-d14784f4aa12
md""" Resuelve numéricamente el siguiente problema de Dirichlet:

$u''(x) + u(x) = 0,\;\;0<x<\pi$

$u(0)=\alpha ,\;\; u(\pi) = \beta$

(a) Si $\alpha = 1$ y $\beta = -1$. ¿A qué solución parece converger cuando 
$h\rightarrow 0$?"""


# ╔═╡ cca6e5b7-fde4-4230-993a-b42482e803f2
md""" # Solución:"""

# ╔═╡ 9899eb31-e9e3-4a40-91e3-81bc220c046f
md""" Uso la función de ecuación general con condiciones de Dirichlet"""

# ╔═╡ bbbb80c4-1ddc-46c6-a651-ada283cbbe64
md"""Para  

$h=\frac{1}{10},\frac{1}{100},\frac{1}{1000}$:"""

# ╔═╡ 9b52689c-cf08-4aca-88a9-6b08363a8fc2
xex3,yex3=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,-1,9);

# ╔═╡ abbc04a0-365a-4360-80ff-296610827f66
xxex3,yyex3=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,-1,99);

# ╔═╡ c3796616-f2e0-45b4-bbc4-e587712c0417
xxxex3,yyyex3=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,-1,999);

# ╔═╡ fab34cb2-2cc8-4dfe-8e61-1fa9800a6ad2
plot([xxxex3, xxex3, xex3],[yyyex3, yyex3, yex3],title = "β = -1",label = ["h = 0.1" "h = 0.01" "h = 0.001"])

# ╔═╡ 1716347f-297e-4a06-b3bb-486b688d7e65
md""" (b) Cambiando a $\beta = 1$, ¿cómo se comporta la solución cuando $h\rightarrow 0$?"""

# ╔═╡ 65c7fe4e-b0a6-4bc4-b4a6-e01dcf93f25e
md""" # Solución:"""

# ╔═╡ 3ebe6eed-802d-4c27-b3c3-442e084438cc
xex3b,yex3b=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,1,9);

# ╔═╡ f6c6bfab-99e1-4782-bbb4-e0d03d75b8ae
xxex3b,yyex3b=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,1,99);

# ╔═╡ 678f0bce-ae41-4411-b3c3-7bce7150796e
xxxex3b,yyyex3b=EGCD(x->1,x->0,x->1,x->0,[0,pi],1,1,999);

# ╔═╡ 03128ee5-a47c-4b9d-8291-2c775a391832
plot([xxxex3b, xxex3b, xex3b],[yyyex3b, yyex3b, yex3b],title = "β = 1",label = ["h = 0.1" "h = 0.01" "h = 0.001"])

# ╔═╡ 5d5eae91-3191-4778-8f01-d93f81f5c506
md""" # Ejercicio 4:"""

# ╔═╡ 5e2214f5-2b65-47fc-a272-ab0ba0341909
md""" Considera la solución numérica del siguiente problema de perturbaciones singulares:

$\epsilon u''(x) - u'(x) = f(x)$

$u(0)=\alpha ,\;\; u(1) = \beta$

Donde $\epsilon>0$ es un parámetro pequeño llamado el número de Péclet. Dibuja las soluciones para $\alpha=1$, $\beta=3$ y $f(x)=-1$, para distintos valores de $\epsilon=0.3,0.1,0.05,0.01$."""


# ╔═╡ ba53ca83-70ae-4c70-8e28-d7ce3b8e76e2
md""" # Solución:"""

# ╔═╡ 5cff3ae9-81ee-4acd-9663-6b9b8ed6782d
md""" Tenemos que $\lambda(x) = \epsilon$, $\mu(x) = -1$ y $\eta(x) = 0$ luego con la función de función general con condiciones de Dirichlet tenemos y tomando $m=100$:"""

# ╔═╡ 907fcbea-6e46-4609-abe8-af42f7adeb36
x4,y4 = EGCD(x->0.3,x->-1,x->0,x->-1,[0,1],1,3,100);

# ╔═╡ c49c8c5d-4c45-42e6-9e23-a958445ef6f6
xx4,yy4 = EGCD(x->0.1 ,x->-1,x->0,x->-1,[0,1],1,3,100);

# ╔═╡ f7d443c8-cfe2-4fa7-903c-9a08a552bdea
xxx4,yyy4 = EGCD(x->0.05 ,x->-1,x->0,x->-1,[0,1],1,3,100);

# ╔═╡ e47a0593-6ecd-4944-bf8c-5eb8dca295a7
xxxx4,yyyy4 = EGCD(x->0.01 ,x->-1,x->0,x->-1,[0,1],1,3,100);

# ╔═╡ 4692d528-e99b-4ba9-a2fc-f097821620fe
plot([xxxx4, xxx4, xx4, x4],[yyyy4, yyy4, yy4, y4],label = ["ϵ = 0.01" "ϵ = 0.05" "ϵ = 0.1" "ϵ = 0.3"])

# ╔═╡ ea52973e-899b-425e-8707-b5e741099fe3
md""" # Ejercicio 5:"""

# ╔═╡ 933f6417-94a7-4b16-8aac-5b6cfb22759e
md""" Considera el problema de contorno

$u''(x)=e^{x},\;\;0\leq x\leq 3,$

$u(0)=-5,\;\;u(3)=3$

(a) Determina el conjunto de nodos de Chebyshev definidos en la teoría y guarda los valores en una variable xnod."""

# ╔═╡ 78b4541c-83aa-45cb-8ffc-e1e50725acaa
md""" # Solución:"""

# ╔═╡ f017bfb8-de47-4580-b340-aed490c56444
md""" Defino función que me crea nodos de Chebyshev"""

# ╔═╡ 1742702b-5eaa-4bb5-a2ff-d26388594aee
function Cheb_nod(intervalo, m)
	z=zeros(m+2)
	x=zeros(m+2)
	for i in 0:(m+1)
		z[i+1] = i/(m+1)
		x[i+1] = intervalo[1] + (1/2)*(intervalo[2]-intervalo[1])*(1+cos(pi*(1-z[i+1])))
	end
	return x
end

# ╔═╡ 5639133b-41a7-45fd-8038-569f2a9e18e8
xnod=Cheb_nod([0,3],5)

# ╔═╡ 8220597c-9538-49c4-9d57-9adac3cc37a3
md""" (b) Utilizando la función fdcoeff definida en la Práctica 1, calcula la matriz 
 del sistema de diferencias finitas para el conjunto de nodos xnod."""

# ╔═╡ 5a921555-cba6-4f4b-bedb-c584d05f675d
md""" # Solución:"""

# ╔═╡ e50705d7-d02e-4ed7-8a4f-999c674f38b6
md""" Cargo la función fdcoeff"""

# ╔═╡ 7bb9b6c3-839e-48b8-882e-08dc18271732
function fdcoeff(k,x,nodos)
	nnodos = length(nodos)
	A = ones(nnodos,nnodos)
	Δx = nodos .-x
	for i = 2:nnodos
		A[i,:] = Δx.^(i-1)/factorial(i-1)
	end
	y = zeros(nnodos)
	y[k+1] = 1
	return A\y #así se resuelve el sistema
end

# ╔═╡ 903a432c-67a0-4ce3-8ba4-6dec27f63722
md""" Voy a crear una función para resolver el sistema $Ax=b$ con nodos de Chebyshev"""

# ╔═╡ 8a823b02-488a-4bc8-850e-a24f066b7121
function Cheb_met(f,a1,b1,alpha,beta,m)

	b = zeros(m+2)
	z = zeros(m+2)
	x = zeros(m+2)
	A = zeros(m+2,m+2)
	
	for k =0:m+1
		z[k+1] = k/(m+1)
		x[k+1] = a1 + (1/2)*(b1-a1)*(1 + cos(pi*(1-z[k+1])))
	end
	
	b[1] = alpha
	b[m+2] = beta
	
	A[1,:] = fdcoeff(0,x[1],x)
	for i = 2:m+1
		b[i] = f(x[i])
		A[i,:] = fdcoeff(2,x[i],x)
	end
	A[m+2,:] = fdcoeff(0,x[m+2],x)
	
	return x, A\b
end

# ╔═╡ a529b6ae-0e44-46ee-ae3d-b898c258fb88
md""" c) Resuelve el sistema para distintos valores de $m=5,7,9,11,13$ y haz una gráfica para el error del máximo (norma infinito) versus $h=1/(m+1)$ a escala logarítmica, teniendo en cuenta que la solución exacta viene dada por

$e^{x} + (9 - e^{3})\frac{x}{3} - 6$

Comprueba que la pendiente no es constante y es cada vez más alta a medida que $h$ se hace más pequeño."""

# ╔═╡ d9f76cde-8e4b-4ad9-be05-43a41661ce00
xex5,yex5 = Cheb_met(x->exp(x),0,3,-5,3,5);

# ╔═╡ 7c9c866e-f2f2-4bf8-a1a7-bd8312e7be92
xex55,yex55 = Cheb_met(x->exp(x),0,3,-5,3,7);

# ╔═╡ 95c860be-0937-4d4f-a05d-e5d57e575788
xex555,yex555 = Cheb_met(x->exp(x),0,3,-5,3,9);

# ╔═╡ 2296f778-8802-48ca-bdeb-bceddf602a45
xex5555,yex5555 = Cheb_met(x->exp(x),0,3,-5,3,11);

# ╔═╡ e8ee5dd7-9c8e-4abe-aa23-ea4ffbfbc8ea
xex55555,yex55555 = Cheb_met(x->exp(x),0,3,-5,3,11);

# ╔═╡ 67eecae0-8e4a-40a9-a04b-20d6515ecd5e
function sole(x)
	sol = exp(x) + (9 - exp(3))*(x/3) - 6
	return sol
end

# ╔═╡ 04c23685-d6bf-4634-a29a-af9c68dd56ca
begin
	ss=zeros(length(xex55555));
	ss=sole.(xex55555);
end

# ╔═╡ 34153f67-343d-4fa8-ab38-31606be6b8af
plot(xex55555,[yex55555, ss], label = ["Solución numérica" "Solución exacta"], title = "Método espectral de Chebyshev")

# ╔═╡ 6283e263-48d5-433a-9488-488d90ebf992
md""" Calculemos el error"""

# ╔═╡ 874d31f9-d150-455e-a130-61a7f795c7fc
begin
	e=zeros(5);
	h=zeros(5);
end

# ╔═╡ fb81e8ba-9d9c-4a00-8bb4-3a409831bc3f
for i in 5:2:13
	t=Int((i-5)/2)
	h[t+1]=1/(i+1)
	diff= Cheb_met(x->exp(x),0,3,-5,3,i)[2] -        
	sole.(Cheb_met(x->exp(x),0,3,-5,3,i)[1])
	e[t+1]=maximum(abs.(diff))
end

# ╔═╡ 8674606f-9c47-4810-9a1f-e6c0b0d3f026
plot(h,e,xaxis =:log,yaxis=:log, xlab="log h",ylab = "log errmax", title = "Método espectral de Chebyshev")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
Plots = "~1.38.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "619af3f8b197e2fa7d32409e4f6493e462127418"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "660b2ea2ec2b010bb02823c6d0ff6afd9bdc5c16"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d5e1fd17ac7f3aa4c5287a61ee28d4f8b8e98873"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "071602a0be5af779066df0d7ef4e14945a010818"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.22"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6f4fbcd1ad45905a5dee3f4256fabb49aa2110c6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.7"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "8ac949bd0ebc46a44afb1fdca1094554a84b086e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.5"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6edfe154ad7b313c01aceca188c05c835c67360"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.4+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═6892b420-b0fe-11ed-1458-711f19c46668
# ╠═159bb200-969e-4928-a6d8-9e43fe7d07ac
# ╠═cf2dfedc-f412-4261-a5d9-9e557eed5c42
# ╠═eb9c4ea1-308c-424c-8b1d-7b936f68b66e
# ╠═a39b572a-3e38-4b10-b215-1855555b1cd4
# ╠═d4de887a-615a-4d38-8b62-5946b27a9185
# ╠═743a19df-7b68-48bb-abe2-1480a445e7e8
# ╠═266b104c-c0ec-4fa5-b706-7d2566f3c531
# ╠═28dd18e2-6651-4e32-9af0-682c5ce9d521
# ╠═d5345ec5-2771-49a6-83dc-79e598493460
# ╠═e0050aa5-5403-418f-b44e-e9bd88b51f88
# ╠═b5520aa5-63ec-41ff-9f4f-20ddd85d5e9c
# ╠═a206d82e-010c-486a-87be-f90bd4c114ef
# ╠═10d5d1e2-e085-4acb-8e90-50799cb91fe0
# ╠═4a7679e7-edf6-4ad3-8397-326c05f683e4
# ╠═1cd9579e-59b8-49ac-886b-3b24b1faf030
# ╠═49ff9dbe-e0da-4072-8dcd-c693feaad5f0
# ╠═42724883-c2b4-4d3b-b2c3-b9819be1b5aa
# ╠═f61fe6d6-7724-4368-bcc7-1c28269318f3
# ╠═5690b3cb-8fbc-46f3-9a4e-e8e7be1b8919
# ╠═f29616f9-459f-4ddc-803c-04eea8d41f48
# ╠═fad3e8ba-5797-42a4-aa13-59d7f5033eb0
# ╠═ac812159-7bc1-4fea-bd1d-323edc85b1cf
# ╠═ba23629d-0664-4af3-aba0-cebcc2e332a9
# ╠═41246027-3220-4256-bd2e-f8c5f34c8176
# ╠═bdc568c3-16ee-4b9f-ad42-a53067e8ebe2
# ╠═b2b6d7dc-b9d1-4b94-8af0-624d59e72e17
# ╠═704af3a8-5836-4400-8cd9-6b30efa514ee
# ╠═35faa2fd-e37c-4f9a-87cb-90fe2513d049
# ╠═9b01627c-b6a4-4865-96b0-9f6483616eb4
# ╠═5b7e8797-2b96-4ee7-ae72-4a333369d289
# ╠═8d9dd22f-d8cd-4720-8248-e4b58dd47a83
# ╠═acd29996-cfc4-4e16-baae-5ffc6ea26a47
# ╠═7bd417ef-db85-472d-bcd9-0f2a5f8f43ba
# ╠═70f3183b-2818-4249-8414-6a396f0d073f
# ╠═ed4365e0-205a-4f9d-b2a0-cb22ba65c730
# ╠═3503d295-7f6b-4f9b-848b-6159317fd19a
# ╠═c1741f3d-9e66-4b10-9e37-d3606842f704
# ╠═6b87c10e-1856-47ff-994a-d83bc6d01ace
# ╠═67c329be-ea7f-4dde-bd80-7fe556f7a74a
# ╠═fd007c48-34fc-4fee-b936-f66ce9f6128e
# ╠═f9ad13d1-0645-4482-813a-dfdcd16ac247
# ╠═43b1bb5f-a431-4a7f-980c-150ea0db6473
# ╠═cfa02242-2d07-4e40-aea5-d2dd8f5032dd
# ╠═ebc4c77e-025c-416a-9b0d-f4e7de796928
# ╠═9de8ba60-db2f-41a7-a31d-7579f5158bd0
# ╠═23e684b0-f03d-4a92-aee2-4d63960ba796
# ╠═8d938768-d418-4b89-a9b7-d27c9f748730
# ╠═8dc14c5e-b6a3-480a-9d30-3476805557ee
# ╠═b69d1fe6-445c-4fb8-9a4c-a245799b3e28
# ╠═cb05c42b-afb0-4c73-aac1-8ebfcacc064b
# ╠═e1810634-1ff5-43c4-bf81-5e6870c07a25
# ╠═3da56891-d0ce-4287-9e80-4105a6da99ff
# ╠═379238de-646c-4db9-a853-7593616bd1c9
# ╠═4ad30e26-755d-4cf1-a5f3-af714580ab5f
# ╠═fe5cb681-30c8-4dca-8744-955a80cb1058
# ╠═ffcb9628-d626-42c2-8519-d9df668663df
# ╠═e9068806-df9f-4990-987a-792a3fa9d34d
# ╠═f2f75266-6124-4656-bdbe-8bdaecb02cf0
# ╠═5240a2a5-2da0-4c66-9d48-c9baf3762e29
# ╠═b8376f9f-a95a-4bb8-bb00-fd149649be91
# ╠═03539f0a-b958-49d9-a80e-8c94ebae688e
# ╠═da37f334-8b65-46cd-b085-9dca5bfe6d3f
# ╠═7859444f-9afb-4e24-8909-e6a2665e5242
# ╠═e3e5ff04-8683-4410-8204-bda2807a9e8b
# ╠═46bbccc1-1a57-419c-b30d-30e24bea8cf4
# ╠═6bed92ed-1234-4081-b386-29cc7eb1281a
# ╠═3e4c1677-bee6-485e-bd20-08d68bb7972f
# ╠═02f1147a-7995-429b-a56d-564548e649f9
# ╠═7319f969-ae44-4363-8b65-d14784f4aa12
# ╠═cca6e5b7-fde4-4230-993a-b42482e803f2
# ╠═9899eb31-e9e3-4a40-91e3-81bc220c046f
# ╠═bbbb80c4-1ddc-46c6-a651-ada283cbbe64
# ╠═9b52689c-cf08-4aca-88a9-6b08363a8fc2
# ╠═abbc04a0-365a-4360-80ff-296610827f66
# ╠═c3796616-f2e0-45b4-bbc4-e587712c0417
# ╠═fab34cb2-2cc8-4dfe-8e61-1fa9800a6ad2
# ╠═1716347f-297e-4a06-b3bb-486b688d7e65
# ╠═65c7fe4e-b0a6-4bc4-b4a6-e01dcf93f25e
# ╠═3ebe6eed-802d-4c27-b3c3-442e084438cc
# ╠═f6c6bfab-99e1-4782-bbb4-e0d03d75b8ae
# ╠═678f0bce-ae41-4411-b3c3-7bce7150796e
# ╠═03128ee5-a47c-4b9d-8291-2c775a391832
# ╠═5d5eae91-3191-4778-8f01-d93f81f5c506
# ╠═5e2214f5-2b65-47fc-a272-ab0ba0341909
# ╠═ba53ca83-70ae-4c70-8e28-d7ce3b8e76e2
# ╠═5cff3ae9-81ee-4acd-9663-6b9b8ed6782d
# ╠═907fcbea-6e46-4609-abe8-af42f7adeb36
# ╠═c49c8c5d-4c45-42e6-9e23-a958445ef6f6
# ╠═f7d443c8-cfe2-4fa7-903c-9a08a552bdea
# ╠═e47a0593-6ecd-4944-bf8c-5eb8dca295a7
# ╠═4692d528-e99b-4ba9-a2fc-f097821620fe
# ╠═ea52973e-899b-425e-8707-b5e741099fe3
# ╠═933f6417-94a7-4b16-8aac-5b6cfb22759e
# ╠═78b4541c-83aa-45cb-8ffc-e1e50725acaa
# ╠═f017bfb8-de47-4580-b340-aed490c56444
# ╠═1742702b-5eaa-4bb5-a2ff-d26388594aee
# ╠═5639133b-41a7-45fd-8038-569f2a9e18e8
# ╠═8220597c-9538-49c4-9d57-9adac3cc37a3
# ╠═5a921555-cba6-4f4b-bedb-c584d05f675d
# ╠═e50705d7-d02e-4ed7-8a4f-999c674f38b6
# ╠═7bb9b6c3-839e-48b8-882e-08dc18271732
# ╠═903a432c-67a0-4ce3-8ba4-6dec27f63722
# ╠═8a823b02-488a-4bc8-850e-a24f066b7121
# ╠═a529b6ae-0e44-46ee-ae3d-b898c258fb88
# ╠═d9f76cde-8e4b-4ad9-be05-43a41661ce00
# ╠═7c9c866e-f2f2-4bf8-a1a7-bd8312e7be92
# ╠═95c860be-0937-4d4f-a05d-e5d57e575788
# ╠═2296f778-8802-48ca-bdeb-bceddf602a45
# ╠═e8ee5dd7-9c8e-4abe-aa23-ea4ffbfbc8ea
# ╠═67eecae0-8e4a-40a9-a04b-20d6515ecd5e
# ╠═04c23685-d6bf-4634-a29a-af9c68dd56ca
# ╠═34153f67-343d-4fa8-ab38-31606be6b8af
# ╠═6283e263-48d5-433a-9488-488d90ebf992
# ╠═874d31f9-d150-455e-a130-61a7f795c7fc
# ╠═fb81e8ba-9d9c-4a00-8bb4-3a409831bc3f
# ╠═8674606f-9c47-4810-9a1f-e6c0b0d3f026
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
