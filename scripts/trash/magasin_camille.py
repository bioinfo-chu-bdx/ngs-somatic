#!/usr/bin/env python
# -*- coding: utf-8 -*-

def afficher_stock(stock):
	print "\nArticles en stock :"
	for article in stock:
		print "- %s (prix : %s; quantité : %s)" % (article, stock[article]["prix"], stock[article]["quantite"])

def selection_article(stock): # fonction permettant de selectionner un article a vendre, et verifie s'il existe s'il est toujours en stock
	afficher_stock(stock)
	selection = raw_input("\n> Tapez le nom de l'article à encaisser :\n")
	if not selection in stock :
		print "\nErreur: l'article %s n'existe pas" % (selection)
	elif stock[selection]["quantite"] <= 0:
		print "\nErreur: l'article %s n'est plus en stock" % (selection)
	else:
		return selection

def afficher_caisse(caisse):
	print "\nEtat de la caisse :"
	for espece in sorted(caisse.keys(),reverse=True): # trie du plus grand (billet de 50) au plut petit (piece de 1 centime)
		if espece >= 5:
			print "- %s billets de %s euros" % (caisse[espece], espece) # si au dessus de 5, on parle de billets
		else:
			print "- %s pièces de %s euros" % (caisse[espece], espece) # si en dessous de 5, on parle de pieces
	total = calcul_total_caisse(caisse)
	print "\nTotal dans la caisse : %s euros" % (total)
			
def calcul_total_caisse(caisse):
	total = 0
	for espece in caisse:
		total = total + (espece * caisse[espece])
	return total

def encaisser(caisse,prix_article): # fonction permettant de calculer les billets et pieces à rendre
	encaissement_reussi = False
	total_caisse = calcul_total_caisse(caisse)
	montant_recu = False
	while not montant_recu:
		try:
			recu = raw_input("\n> Tapez le montant donner par le client :\n")
			recu = recu.replace(',','.') # si l'utilisateur tape une virgule au lieu d'un point
			recu = float(recu)
		except:
			print "Erreur: format non-valable. Saisir le montant en numérique (ex : 20.50)."
			continue
		if recu < prix_article:
			print "Erreur: Montant inférieur au prix de l'article."
		else:
			montant_recu = recu
			
	a_rendre = montant_recu - prix_article
	if a_rendre > total_caisse:
		print "Erreur : pas assez d'argent de la caisse pour rendre la monnaie au client"
	else:
		# Encaissement du montant recu (on estime le detail en billet et pieces)
		print "\nRajout dans la caisse de :"
		for espece in sorted(caisse.keys(),reverse=True):
			while ((montant_recu > espece) or (str(montant_recu) == str(espece))):
				caisse[espece] = caisse[espece] + 1
				print "\t- %s euros" % espece
				montant_recu = montant_recu - espece
		# Rendu de la monnaie
		print "Rendu au client de %s euros :" % a_rendre
		for espece in sorted(caisse.keys(),reverse=True):
			while ((a_rendre > espece) or (str(a_rendre) == str(espece))):
				caisse[espece] = caisse[espece] - 1
				print "\t- %s euros" % espece
				a_rendre = a_rendre - espece
		encaissement_reussi = True
	return encaissement_reussi
				
#########################################################
				
caisse = {
	500.0:0,
	200.0:0,
	100.0:0,
	50.0:4,
	20.0:6,
	10.0:10,
	5.0:10,
	2.0:10,
	1.0:20,
	0.50:15,
	0.20:30,
	0.10:30,
	0.05:20,
	0.02:30,
	0.01:30
}

stock = {
	"Chemise noire" : {"prix" : 40.0, "quantite" : 6},
	"Chemise bleue" : {"prix" : 40.0, "quantite" : 2},
	"Pantalon kaki" : {"prix" : 59.90, "quantite" : 2},
	"Jean bleue" : {"prix" : 60.50, "quantite" : 5},
	"Casquette" : {"prix" : 19.90, "quantite" : 3},
	"Baskets adidas Pharrell Williams" : {"prix" : 41.99, "quantite" : 1},
}

print "Caisse du magasin bonjour."

fermer_caisse = False
while not fermer_caisse:
	action = raw_input("\n> Selectionner l'action à effectuer : \n\t 1-Vendre un article \n\t 2-Afficher stock \n\t 3-Afficher caisse \n\t 4-Quitter\n")
	
	if action == "1":
		# Selection de l'article
		article = False
		while article == False:
			article = selection_article(stock)
			if not article :
				print "\nRetour à la selection de l'article."
				article = False
			else:
				print "\nArticle selectionné : %s" % (article)

		# Encaissement
		prix_article = stock[article]["prix"]
		print "prix de l'article : %s euros" % prix_article
		encaissement = encaisser(caisse,prix_article)
		
		if encaissement:
			# Article vendu, quantité stock -1
			stock[article]["quantite"] = stock[article]["quantite"] -1
		else:
			print "Erreur lors de l'encaissement. Annulation de la transaction." 
			print "Retour au menu principal."
		
	elif action == "2":
		afficher_stock(stock)
	elif action == "3":
		afficher_caisse(caisse)
	elif action == "4":
		fermer_caisse = True
	else:
		print "\nErreur: s.v.p. tapez 1, 2, 3 ou 4\n"

print "Bye-Bye"
exit()
