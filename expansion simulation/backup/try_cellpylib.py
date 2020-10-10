import cellpylib as cpl
# Glider
cellular_automaton = cpl.init_simple2d(60, 60)
cellular_automaton[:, [28,29,30,30], [30,31,29,31]] = 1
# Blinker
cellular_automaton[:, [40,40,40], [15,16,17]] = 1
# Light Weight Space Ship (LWSS)
cellular_automaton[:, [18,18,19,20,21,21,21,21,20],
[45,48,44,44,44,45,46,47,48]] = 1
# evolve the cellular automaton for 60 time steps
cellular_automaton = cpl.evolve2d(cellular_automaton, timesteps=60,
neighbourhood='Moore',
apply_rule=cpl.game_of_life_rule)
cpl.plot2d_animate(cellular_automaton)

#%% rule table
import cellpylib as cpl
rule_table, actual_lambda, quiescent_state = cpl.random_rule_table(lambda_val=0.45, k=4, r=2,
                                                                   strong_quiescence=True,
                                                                   isotropic=True)
cellular_automaton = cpl.init_random(128, k=4)
# use the built-in table_rule to use the generated rule table
cellular_automaton = cpl.evolve(cellular_automaton, timesteps=200,
                                apply_rule=lambda n, c, t: cpl.table_rule(n, rule_table), r=2)
