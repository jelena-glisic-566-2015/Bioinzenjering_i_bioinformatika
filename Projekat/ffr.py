import math
alfa_tabela = [1,0.9,0.8, 0.6,0.4,0.2,0.1,0.05]
coeff_tabela = [0, 0.0002, 0.0008, 0.0034, 0.0087, 0.0182, 0.0264, 0.0348]

#Simpsonovo pravilo aproksimacije integrala
#f - funkcija, a,b - interval, n - broj segmenata

def simpson(f, a, b, n):
    h=(b-a)/n
    k=0.0
    x=a + h
    for i in range(1,n/2 + 1):
        k += 4*f(x)
        x += 2*h

    x = a + 2*h
    for i in range(1,n/2):
        k += 2*f(x)
        x += 2*h
    return (h/3)*(f(a)+f(b)+k)

def interpolacija(coeff):
	if coeff > 0.0348: return 0.049 #vrati nesto manje od 0.05, da bi moglo u stablu odlucivanja da se uzme prava jednacina
	if coeff <= 0 : return 1
	#pronadji najblizu vrednost koeficijenta iz tabele
	coeff_tab = min(coeff_tabela, key = lambda k: abs(k-coeff))
	najblize_alfa = alfa_tabela[coeff_tabela.index(coeff_tab)]
	#da li uzeti sledecu ili proslu vrednost za interpolaciju 
	direction = coeff - coeff_tab
	if direction > 0:
		try: 
			index = coeff_tabela.index(coeff_tab)+1
			key = coeff_tabela[index]
			value = alfa_tabela[index] 
			return najblize_alfa + (coeff - coeff_tab) * (value - najblize_alfa)/(key-coeff_tab)
		except IndexError: #poslednji clan
			index = coeff_tabela.index(coeff_tab)-1
			key = coeff_tabela[index]
			value = alfa_tabela[index]
			return value + (coeff - key) * (najblize_alfa - value)/(coeff_tab-key)
	elif direction < 0:
		try: 
			index = coeff_tabela.index(coeff_tab)-1
			key = coeff_tabela[index]
			value = alfa_tabela[index]
			return value + (coeff - key) * (najblize_alfa - value)/(coeff_tab-key)
		except IndexError: #prvi clan
			index = coeff_tabela.index(coeff_tab)+1
			key = coeff_tabela[index]
			value = alfa_tabela[index]
			return najblize_alfa + (coeff - coeff_tab) * (value - najblize_alfa)/(key-coeff_tab)

def func(x):
	return float(1+4*x+9*x*x+4*x*x*x)/float(x*(3+2*x)*(3+2*x+x*x)*(3+2*x+x*x)) #A3 jednacina

def func1(x):
	return (1-x)/5 * func(x) * (6+x) #A2 jednacina

def func2(x):
	return 8*math.pi*0.00278*Q/(CSAst*CSAst) #A1 jednacina

l_stenosis = float(raw_input("Uneti duzinu stenoze [mm]: "))
L_sud = float(raw_input("Uneti duzinu suda [mm]: "))
Q = float(raw_input("Uneti protok krvi [ml/s]: "))
CSAst = float(raw_input("Uneti poprecni presek stenoze [mm^2]: "))
CSAdist = float(raw_input("Uneti distalni poprecni presek [mm^2]: "))
CSAin = float(raw_input("Uneti poprecni presek ulaza [mm^2]: "))
CSAout = float(raw_input("Uneti poprecni presek izlaza [mm^2]: "))
Rst = math.sqrt(CSAst)/math.pi
Pa = float(raw_input("Pritisak aorte [mmHg]: "))#pretvaranje mmHg u Pa, pritisak sredji arterijski
kinematska_viskoznost = 2.65
rho = 10.60
#koeficijent L/ReD
coeff = math.pi * kinematska_viskoznost * l_stenosis / (4*Q*1000) #Q*1000 jer pretvaramo ml u mm^3
print "Koeficijent = " + str(coeff)
alfa = interpolacija(coeff)
print "alfa = " + str(alfa) 
	
deltaP = (float(rho*Q*Q)/(2*CSAst*CSAst)) * 96.0/5. * simpson(func, alfa, 1, 1000) #malo p

print "FFR: " + str(100 * (Pa - deltaP)/Pa) 
Pstenosis = 0
Puniform = rho*Q*Q/2*(1/CSAst - 1/CSAdist)*(1/CSAst - 1/CSAdist)
Pparabolic = rho*Q*Q*(1/CSAst - 1/CSAdist)*(1/CSAst - 1.0/3 * 1/CSAdist)
if alfa < 0.05:
	L_ulaz = Q/(math.pi*kinematska_viskoznost)*simpson(func1, alfa, 1, 1000)
	Pstenosis = float(rho*Q*Q)/(2*CSAst*CSAst) * 96/5 * simpson(func, alfa, 1, 1000) + simpson(func2, 0, L_ulaz-l_stenosis, 1000) + Pparabolic
	print "Pad pritiska kroz stenozu: " + str(Pstenosis) 
	DeltaP = rho*Q*Q/2* (1./(CSAout*CSAout) - 1./(CSAin*CSAin))+ float(rho*Q*Q)/(2*CSAst*CSAst) * 96/5 * simpson(func, alfa, 1, 1000) +  simpson(func2, 0, L_ulaz-l_stenosis, 1000)  + rho*Q*Q*(1./CSAst - 1./CSAdist )*(1/CSAst - 1.0/3 * 1/CSAdist) #veliko P
	print "Ukupni pad pritiska: " + str(DeltaP)
else:
	Pblunt = Puniform + (Pparabolic - Puniform)*(1-alfa)*(1-alfa)
	Pstenosis = float(rho*Q*Q)/(2*CSAst*CSAst) * 96/5 * simpson(func, alfa, 1, 1000)+Pblunt
	print "Pad pritiska kroz stenozu: " + str(Pstenosis) 
	DeltaP =  rho*Q*Q/2* (1./(CSAout*CSAout) - 1./(CSAin*CSAin))+ float(rho*Q*Q)/(2*CSAst*CSAst) * 96/5 * simpson(func, alfa, 1, 1000) +  simpson(func2, 0, L_sud-l_stenosis, 1000) + rho*Q*Q/2* ((1/CSAst - 1/CSAdist)*(1/CSAst - 1/CSAdist) + ( 2* (1./CSAst - 1./CSAdist )*(1/CSAst - 1.0/3 * 1/CSAdist) - (1/CSAst - 1/CSAdist) * (1/CSAst - 1/CSAdist)) * (1-alfa) * (1-alfa)) #veliko P
	print "Ukupni pad pritiska: " + str(DeltaP)
	
