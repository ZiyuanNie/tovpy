# Equation of state tables

EOS tables in format 

    L1
    C1 C2 C3 C4
    ... 

with 

    L1 number of tabulated points
    C1 energy density/c^2    (g/cm^3)
    C2 pressure              (dynes/cm^2)
    C3 enthalpy              (cm^2/s^2)
    C4 baryon number density (cm^{-3})

This format matches what is used in Stergioulas' [rns code](https://github.com/cgca/rns).
Some EOS tables are taken from there:

     eosA		PANDHARIPANDE NEUTRON: A&B EOS A
     eosAU		WFF1 (denoted AV14+UVII in WFF) matched to Negele-Vautherin
     eosB		PANDHARIPANDE HYPERON: A&B EOS B
     eosC		BETHE-JOHNSON MODEL 1: A&B EOS C
     eosF		PANDHARIPANDE (HYPERON) AND ARPONEN: A&B EOS F
     eosFP		FRIEDMAN&PANDHARIPANDE EOS.(NEGELE-VAUTHERIN FOR N<.1FM**-3,THEN * BPS FOR N<.001FM**-3)
     eosFPS	        Lorenz, Ravenhall and Pethick, 1993, PRL 70,379
     eosG		CANUTO-CHITRE (SOLID): A&B EOS G 
     eosL		PANDHARIPANDE AND SMITH (MEAN FIELD): A&B EOS L
     eosN		WALECKA-SEROT: A&B EOS N
     eosNV		Negele & Vautherin, 1973, Nucl. Phys. A207, 298
     eosO		BOWERS, GLEESON, AND PEDIGO. A&B EOS O
     eosUU		WFF2 (denoted UV14+UVII in WFF) matched to Negele-Vautherin
     eosWNV		WFF3 (denoted UV14+TNI in WFF) matched to Negele-Vautherin
     eosWS		WFF3 (UV14+TNI) mathed to EOS FPS

A&B = Arnett and Bowers (1977) APJS 33, 415
WFF = Wiringa, Fiks and Fabrocini (1988), Phys. Rev. C 38, 1010

