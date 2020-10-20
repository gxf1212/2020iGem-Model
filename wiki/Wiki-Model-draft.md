---
author:
- Xufan Gao
bibliography:
- './iGem.bib'
title: |
    Wiki-Model
---

# 要求

快去借鉴往年的！！

## 格式

- 每个模型一个网页，故用一级标题
- 现在字太多了。。。。介绍and结果：画漂亮的图，当美工（安排学习）
- 尾部添加附录，下载or链接or展示，源代码
- 假设的首句、biomass concentration等词，给个颜色或字体？我先黑体了
- 

## 内容

- 模块化，然后可随便排版。所以两种结构
  - 背景介绍、主要假设、数学表达&参数、结果分析
  - 背景+假设，数学表达、参数+部分假设（结合，能说清）、结果分析
- 大家也可写，要求在grammarly等上修改润色。对我个人写的，则从读者的角度出发简化
- 突出思想，在别人能听懂的地方详细，给人“我看懂了，你们想的真细”的感觉；我们怎么确定参数，两句话即可，突出严谨性。比如生长模型，突出$f$和资源的完备性，结果突出能共生、结论
- 时态
  - 背景：该啥就啥，可以一现、过去、现完
  - 一般事实（假设、分析图、评价结果等）用一般现在时，we干了什么这些动作（simulated, graphed, asjusted parameters, concluded）用一般过去时，评价还可将来
  - 其他具体看



# Overview



# The Growth Model

## Background and Assumptions

As the literature indicated, bacteria and cyanobacteria are living together in the soil, both known as crucial microorganisms during the development of biological soil crust. Bacteria secrete organic acids and powerful exoenzymes that dissolve inorganic/organic nitrogen and phosphorus for cyanobacteria to use; they also make many other essential compounds like vitamins accessible. On the other side, *Nostoc sp.* and even more cyanobacteria can fix CO2 and N2 from the air, thus offering organic matter and extra nitrogen as fundamental needs for bacteria. 参考文献 **Such a mutualism is the central relationship between them that determines their survival in the soil**. In some sense, they are interdependent in the nutrition-lacking soil. Although some cyanobacteria and *Bacillus* in lakes may create toxins to repress the growth of each other, many studies have proved that *Nostoc sp.* and *Bacillus subtilis* do coexist in the soil and cooperate to make sandy soil a suitable environment for lives. 参考文献

互补营养的示意图or动态

Our growth model focuses on the above relationship, emphasizing the nutritional constraint. For simplicity,  we also make the following principal assumptions:

1. **We ignore other organisms and most environmental changes** (such as rain, wind, $\ce{CO2}$ and or so), aiming to simulate merely the interations between *Bacillus Subtilis* and *Nostoc sp*. Since the soil environment is really complicated, it is impossible to take everything into account. However, in regard of the nitrogen and phosphrus fixing, we may use the average rate considering all microorganisms in the soil.
2. **We assume that the two microorganisms and all related nutrients are distributed evenly in the soil.** In other words, they are living in a solution-like system and spatial heterogeneity is not considered. Thus, we can use biomass concentration and nutrient concentrition to describe this system quantitively.
3. **"Nutrient" means matter that bacteria and cyanobacteria can directly and quickly**, or "utilizable nutrients". C, N, P exist in plenty of forms in the soil. Literature call them "dissolved nitrogen", "soluable phosphorus" and so on. Only low molecular weight (LMW) or ions can be quickly utilized. For example, amino acids and soluable phosphate are "utilizable nutrients", while humus and $\ce{Ca3(PO4)2}$ are not, which need enzymatical decomposition.
4. **We agree that after Nostoc sp. die, the cell matter is quickly decomposed by bacteria and return to "utilizable nutrients", while *Bacillus Subtilis* don't do this.** This is based on studies of 参考文献.
5. Another assumptions of established models we used.

See xxx-link for other assumptions we make.

## Mathematical Model

### Equations 

Main variables of interest are:

- $N_i$ ($\mathrm{g/L}$): **biomass concentration** (cell dry weight per unit volume); 
- $R_j$ ($\mathrm{g/L}$): **nutrient concentrition** (mass of the element per unit volume). 

In the following equations built on classic **resource competition model**, subscript $i$ equals to 1 for *Bacillus subtilis* and 2 for *Nostoc sp.* ; subscript $j$ equals to $c$ for carbon, $n$ for nitrogen and $p$ for phosphorus, respectively. 

| Equation                                                     | Descriptions                                                 |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| $\dfrac{\mathrm{d}N_i}{\mathrm{d}t}=Birth_i - Death_i$ 编号1~5 | For both microorganisms, the net growth rate $\dfrac{\mathrm{d}N_i}{\mathrm{d}t}$ is defined by birth rate minus death rate. |
| $Birth_i=f_i(R) \cdot r_i N_i $                              | Birth rate $Birth_i$ is modified by a **nutrient-constraint factor** $f_i(R)$ which is not more than 1. |
| $Death_i=m_i N_i$                                            | Growth and death rates are proportional to current population $N_i$ quantity with a ratio. |
| $f_i(R)=\min\limits_j\{\dfrac{R_j}{R_j+K_{ij}N_i}\}$         | The constraint factor is built according to **Monod growth kinetics** with Liebig's law of the minimum. |
| $\dfrac{\mathrm{d}R_j}{\mathrm{d}t}=Income_j-\sum\limits_i Q_{ij} Grow_i$ or $\dfrac{\mathrm{d}N_i}{\mathrm{d}t}$ | The concentration of nutrients are dynamically regulated by production and consumption by microorganisms. |
|  6    |  $\dfrac{\mathrm{d}N_i}{\mathrm{d}t}=Birth_i(1-\gamma_i N_1N_2) - Death_i$ | In the later stage of the BSC development, the two microorganisms are abundant to inhibit each other's growth by secreting toxin. We use a model proposed by (Effect of toxic substances on a two-species competitive system). |


Here we will detailedly explain our innovation in the last two equations.

#### nutrient-constraint factor

- $f_i(R)$ looks quite similar to Michaelis-Menten equation, which is used by many former iGem teams. When $R_j$ goes to infinity, $f_i(R)$ approaches to 1, thus have no effect on the birth rate; otherwise, the birth rate seriously decrease, even leading that net growth rate is less than 0. 

- $f_i(R)$ reflects how much each nutrient limits growth. Inspired by Liebig's law of the minimum, the overall effect is determined by the most lacking nutrient. Even if carbon and nitrogen is sufficient, when phosphrous is deficient, microorganisms cannot grow.

- $K_{ij}$ is a half-saturation constant like Michaelis-Menten constant. However, most literatures ignored that the amount of nutrient needed obviously depends on biomass concentration. Therefore, we adopt half-saturation constant **per unit biomass concentration**, forming a variable and more realistic half-saturation constant. This is actually based on 's research.

  > Kinetics of Bacterial Growth: Relationship between Population Density and Specific Growth Rate of Continuous Cultures

#### dynamic change of nutrients

- Based on assumption 4, the loss of nutrients is proportional to net growth rate for *Nostoc sp.* and birth rate for *Bacillus Subtilis*.
- $Q_{ij}$ means consumption rate of nutrient $j$ per unit biomass concentration. Aerobic respiration is alway ongoing, so we have to use the strict ratio between C consumption and biomass increase. As for N and P, we use their proportion in cell dry weight because there is no extra loss.
- The detailed income and output is shown in equation 6 and explained in "parameters".

Accordingly, equation 5 is formulated into:
$$
\begin{align*}
\dfrac{\mathrm{d}R_c}{\mathrm{d}t}&=-Re_1\cdot N_1+NPh\cdot N_2-Q_{1c}Grow_1-Q_{2c}\dfrac{\mathrm{d}N_2}{\mathrm{d}t}\\
\dfrac{\mathrm{d}R_n}{\mathrm{d}t}&=N_{sol}\cdot N_1 + N_{fix}\cdot N_2-Q_{1n}Grow_1-Q_{2n}\dfrac{\mathrm{d}N_2}{\mathrm{d}t}\\
\dfrac{\mathrm{d}R_p}{\mathrm{d}t}&=P_{sol}\cdot N_1-Q_{1p}Grow_1-Q_{2p}\dfrac{\mathrm{d}N_2}{\mathrm{d}t}\\
\end{align*}
$$
​	图？	







6.16e-7  固磷量除以24h内平均生物量  解磷菌的培养条件优化研究




## Result Analysis

### the growth



如何解释菌溶氮速率更大？





### relationship of the parameters

猜想：Q和速率之比 vs 生长率/死亡率



Kp1

- 结论：越小，刚开始藻菌N都小，但三五天（最终）都涨上来

藻



### Optimization of Initial Dosage













## Conclusion and Prospects

- We found that *Bacillus Subtilis* and *Nostoc sp*. can successfully so-survive in the soil maily due to their high unit-population nitrogen fixation rate
- 
- Though being a moderately good simulation of the growth of the two microorganisms, our growth model is not flawless due to neglecting various environmental factors and needs perfection.



# The Suicide Model

## Background

For simplicity,  we make the following principal assumptions:

- We ignore
- We assume that



## Mathematical Model

parameters' values

| Parameter | Descriptions | Value | References |
| --------- | ------------ | ----- | ---------- |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |


Equations and Variables

| Equation | Variables | Descriptions |
| -------- | --------- | ------------ |
|          |           |              |
|          |           |              |
|          |           |              |
|          |           |              |
|          |           |              |
|          |           |              |




## Result Analysis



## Conclusion and Prospects

- 



# The Expansion Model

- Observe the shape野生型生物膜能够通过改变形状和增加生长细胞的比例来打破这个瓶颈，从而加速生长。

## Background

For simplicity,  we make the following principal assumptions:

- We assume that the main force driving bacteria and cyanobacteria to move is the osmotic pressure. Areas with high concentration EPS (and other nutrient molecules & ions) uptake water from the soil and causes swelling of biofilm. This takes not only microorganisms but also nutrients 
- We ignore
- We assume that



## Mathematical Model

parameters' values

| Parameter | Descriptions | Value | References |
| --------- | ------------ | ----- | ---------- |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |
|           |              |       |            |


Equations and Variables

| Equation | Variables | Descriptions |
| -------- | --------- | ------------ |
|          |           |              |
|          |           |              |
|          |           |              |



## Result Analysis



## Conclusion and Prospects

- 



# poster（原文，可借鉴，需改动！）

要求：300词左右，描述（思想）一半多+结论（目的）一半少

Our models focus on three aspects:

- the growth of the two microorganisms
- the function of the suicide gene
- the spread of BSC

## the growth model

In consideration of lack in nutrients and competition between organisms, we are applying nutritional constraint to traditional Malthusian growth model:
$$
\dfrac{\mathrm{d}N_i}{\mathrm{d}t}=\min_j\{\dfrac{R_j}{R_j+K_{ij}N_i}\}\cdot r_iN_i-m_iN_i
$$


The concentration of nutrients are dynamically regulated by consumption due to growth, photosynthesis, phosphorus solubilization etc. , expressed as:
$$
\dfrac{\mathrm{d}R_j}{\mathrm{d}t}=Income_j-\sum_i Consume_{ij} Grow_i
$$

> We also use Luedeking-Piret equation to roughly predict the yield of polysaccharides which is thought to be a indicator of BSC:
> $$
> \dfrac{\mathrm{d}P}{\mathrm{d}t}=\alpha\dfrac{\mathrm{d}N}{\mathrm{d}t}+\beta N
> $$
> where $\alpha$ and $\beta$ are parameters

Finally, we found the two microorganisms successfully so-survive in the soil maily due to their high unit-population nitrogen fixation rate....

We determined the optimal initial biomass concentration of the two microorganisms by xx algorithm.....

## the suicide model

MazF will effectively kill *Bacillus Subtilis* later.









 

EPS production

Because EPS genes are not studied clearly, we use a classical relationship named **Luedeking-Piret equation** to describe how the yield of cell product is affected by population quantity:

 are experimentally estimated parameters, which we can obtain from literatures. 

To accurately figure out how population quantity will change, we consider multiple environmental resources as significant limiting factors. The basic formula is built on classic **resource competition model**:



where growth and death rates are proportional to current population quantity:



Growth rate is modified by environmental resources including carbon source, nitrogen source, phosphorus source. Their limitation effect is represented based on **Monod growth kinetics** with **Liebig's law of the minimum**:



where 

 are measures of resources. They are dynamic and we separate the incoming and outcoming path:



This can reveal the relationship that *Nostoc sp.* produces carbon and nitrogen for *Bacillus Subtilis* and the *Bacillus Subtilis* produces phosphorus back.

Effect of other factors including water (which we assume will be supplied adequately), light intensity, pH, and temperature will be formulated into one term     according to current models.       is multiplied to the growth rate of corresponding microbes.     are less than or equal to 1.

We use these ordinary differential equations to simulate how EPS is produced. The parameters and models are all from literatures. The initial values come either from experiments on our soil samples or literatures.

 

MazF function

MazF protein cleaves mRNA and promotes apoptosis, whose amount is regulated by MazE and arabinose. The MazE/F expression system is:

![img](file:///C:/Users/Lenovo/AppData/Local/Temp/msohtmlclip1/01/clip_image027.gif)

To gain insight into how exocellular arabinose activates cell death, we first consider the transporting systems whose production is also regulated by arabinose:

![img](file:///C:/Users/Lenovo/AppData/Local/Temp/msohtmlclip1/01/clip_image029.jpg)

So the transporting rate is:

![img](file:///C:/Users/Lenovo/AppData/Local/Temp/msohtmlclip1/01/clip_image031.jpg)

Both expression and transporting follow Michaelis–Menten or Hill equation. Besides we have to consider background expression, degradation and cell division.

*mazE/F* genes are promoted by AraE, so we can use some equation to show the effect of ![img](file:///C:/Users/Lenovo/AppData/Local/Temp/msohtmlclip1/01/clip_image033.gif) on expression[[GX1\]](#_msocom_1) . As for outcoming path, MazE is degraded faster than MazF so we ignore MazF’s natural degradation; however, MazF concentration decays with cell dividing. And when MazE and MazF binds (reversibly), the mRNA cleavage mechanism is disabled. Not only the MazE/F complex decays but also MazE was degraded in it. So the equations are:

![img](file:///C:/Users/Lenovo/AppData/Local/Temp/msohtmlclip1/01/clip_image035.jpg)

As soon as the concentration of MazF reaches a threshold, they start to function. 

We use the equations to simulate how much arabinose will start the self-killing process effectively. The simulation will be run on SimBiology platform. 



# 建模wiki借鉴

- https://2019.igem.org/Team:Wageningen_UR/Results/Spatial_Spread
  可参考它的网页格式和动图效果。这是我昨天找到的一个扩散模型，也是动态图。然后他们有一个把模型假设那些隐藏起来的一个按钮吧，看看那些东西好不好做？

- http://2016.igem.org/Team:Chalmers_Gothenburg/Model

  写作参考，情况很像。。

- 鼠标换成图案？

废句子

> where $\alpha$ and $\beta$ are experimentally-estimated or literature-based parameters
>
> to roughly estimate how much *Bacillus Subtilis* will help with exopolysaccharides (EPS) production and soil improvement in the early state of biological soil crust (BSC) development; 