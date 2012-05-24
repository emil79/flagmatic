P = load("44i")
x = polygen(QQ)
K = NumberField(x**3 - 2*x**2 + 2*x - Integer(2)/3, 'x', embedding=RDF(0.5))
x = K.gen()
C = BlowupConstruction(GraphFlag("8:122886633447755156781122334455667788"), weights=[x/4,x/4,x/4,x/4,(1-x)/4,(1-x)/4,(1-x)/4,(1-x)/4], field=K, phantom_edge=(1,4))
