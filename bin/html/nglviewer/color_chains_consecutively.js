   /* http://colorbrewer2.org/#type=qualitative&amp;scheme=Paired&amp;n=9
  Commented colors from chir.ag/projects/name-that-color
*/

   var colorBrewer9Colors = [
       "skyblue", 
       '#fb9a99', 
       '#b2df8a', 
       "tan",     
       '#cab2d6', // Lavender Grey
       '#1f78b4', // Matisse - Darker blue
       '#33a02c', //Apple (green)
       '#e31a1c', //Alizarin Crimson
       '#ff7f00', //Flush orange
       '#707070', //Chris found black 000000 to be too black - so backing off with less dark
                    ]
   var chainColors = ["skyblue","green","magenta","salmon","orange","black"]
   chainColors = colorBrewer9Colors
   function color_chains_consecutively(o,chainColorIndex) {
       chain_color_representations = []
        chainColorIndex = 0;
        o.structure.eachChain(function (cp) {
            if (cp.residueCount > 1) { // Idea is to exclude the ngl endless little one-off hetero chains
                repr = o.addRepresentation("cartoon", {sele: ":" + cp.chainname, color: chainColors[chainColorIndex]});
                chain_color_representations.push(repr)
                chainColorIndex = (chainColorIndex + 1) % chainColors.length
            }
        })
        return chain_color_representations
   };


	function color_chains_by_pathprox_class(disease1or2PathProxScores) {
		let maxScore = -1.0
		let minScore =  1.0
		for (const [key,value] of Object.entries(disease1or2PathProxScores)) {
   		if (value < minScore) {
   		minScore = value;
   		}
   		if (value > maxScore) {
   		maxScore = value;
   		}
		};
		let scoreRange = maxScore - minScore;

   var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
       this.atomColor = function (atom) {
           let atomScoreKey = atom.resno + ":" + atom.chainname;
           let pathProxScore = disease1or2PathProxScores[atomScoreKey];

           let final_color = 0;

           if (pathProxScore < 0)
                     {
                     let whiteIntensity=parseInt((1.0-(pathProxScore/minScore)) * 0xFF)
                     let blueIntensity = 0xFF
                     final_color = (whiteIntensity << 16) |  (whiteIntensity << 8) | (blueIntensity)
                     }
         			else
                     {
                     let whiteIntensity=parseInt((1.0-(pathProxScore/maxScore)) * 0xFF)
                     let redIntensity=0xFF
                     final_color = (redIntensity  << 16) |  (whiteIntensity << 8) | (whiteIntensity)
                     }
                 return final_color;
			   }
		})

		return schemeId;
  };


    	function color_chains_by_pathprox_class(disease1or2PathProxScores) {
		let maxScore = -1.0
		let minScore =  1.0
		for (const [key,value] of Object.entries(disease1or2PathProxScores)) {
   		if (value < minScore) {
       		minScore = value;
       		}
   		if (value > maxScore) {
       		maxScore = value;
       		}
		};
		let scoreRange = maxScore - minScore;

       var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
           this.atomColor = function (atom) {
               let atomScoreKey = atom.resno + ":" + atom.chainname;
               let pathProxScore = disease1or2PathProxScores[atomScoreKey];

               let final_color = 0;

               if (pathProxScore < 0)
                         {
                         let whiteIntensity=parseInt((1.0-(pathProxScore/minScore)) * 0xFF)
                         let blueIntensity = 0xFF
                         final_color = (whiteIntensity << 16) |  (whiteIntensity << 8) | (blueIntensity)
                         }
             			else
                         {
                         let whiteIntensity=parseInt((1.0-(pathProxScore/maxScore)) * 0xFF)
                         let redIntensity=0xFF
                         final_color = (redIntensity  << 16) |  (whiteIntensity << 8) | (whiteIntensity)
                         }
                     return final_color;
			   }
		})

		return schemeId;
      };

	function color_chains_by_alpha_fold_metrics(alpha_fold_metrics) {

       var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
           this.atomColor = function (atom) {
               alpha_fold_metric = parseFloat(alpha_fold_metrics[atom.residueIndex])

               let final_color = 0x000000

               // Colors from alphafold webpage. 
               if (alpha_fold_metric > 90.0) //  Very high (pLDDT > 90) 
                   final_color = 0x0053D6
               else if (alpha_fold_metric > 70.0) // Confident (90 > pLDDT > 70) 
                   final_color = 0x65CBF3
               else if (alpha_fold_metric > 50.0) // Low (70 > pLDDT > 50) 
                   final_color = 0xFFDB13
               else
                    final_color=0xFF7D45 // Very low (pLDDT < 50) 

            return final_color;
			   }
		})

		return schemeId;
      };



	function color_chains_by_rate4site_scores(rate4site_scores) {

       var schemeId = NGL.ColormakerRegistry.addScheme(function (params) {
           // rate4site_scores are normalized to have average=0.0, std=1.0
           // We want to color [-2,2], with all outliers bright 
           this.atomColor = function (atom) {
               let atomScoreKey = atom.resno + ":" + atom.chainname;
               let final_color = 0x000000
               if (! (atomScoreKey in rate4site_scores)) // Unlike pathprox scores which cover every position, it is possible we won't have some scores
                   final_color = 0x707070 // Color those yellow for now
               else {
               let rate4site_score = parseFloat(rate4site_scores[atomScoreKey])

               if (rate4site_score < 0)
                         {
                         let redIntensity=0xFF
                         let whiteIntensity=0x0
                         if (rate4site_score > -2.0) // Lighten a bit inside 2nd stddev
                             whiteIntensity=parseInt((1.0-(rate4site_score/(-2.0))) * 0xFF)

                         final_color = (redIntensity  << 16) |  (whiteIntensity << 8) | (whiteIntensity)
                         }
             			else // rate4site_score >= 0
                         {
                         let blueIntensity = 0xFF
                         let whiteIntensity=0x0
                         if (rate4site_score < 2.0) // Lighten a bit inside 2nd stddev
                             whiteIntensity=parseInt((1.0-(rate4site_score/2.0)) * 0xFF)

                         final_color = (whiteIntensity << 16) |  (whiteIntensity << 8) | (blueIntensity)
                         }
						}

                     return final_color;
			   }
		})

		return schemeId;
      };



