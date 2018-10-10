#### This document is currently a work-in-progress QIIME 2 terminology glossary.(Version 2018.8)
***

**Action**

*method* 或 *visualizer* 的通用术语。  
  

**Artifact**

可以用作QIIME *method* 或 *visualizer* 的输入的数据，或者可以作为QIIME *method* 的输出生成的数据。在写入文件时， *artifacts* 通常具有扩展名```.qza```。
  
**Method**

将 *artifacts* 和 *parameters* 的某种组合作为输入并将一个或多个 *artifacts* 作为输出生成的操作。   
这些输出 *artifacts* 随后可以用作其他QIIME 2 *methods* 或 *visualizers* 的输入。  
*Methods* 可以在QIIME分析中产生中间或终端输出。  
  
  
 
**Parameter**  
对 *actions* 的原始（即，非 *artifacts* ）输入。  
例如，字符串，整数和布尔值是原始输入。  
永远不会从 *action* 输出原生类型编码。  

**Pipeline**  
*Actions* 的组合。  

**Plugin**  
*Plugin* 提供微生物组（特定域的）分析功能，用户可以通过围绕QIIME 2框架构建的各种接口访问这些功能。  
*Plugin* 可以由任何人开发和分发。
在更多技术术语中，插件是Python3包，它实例化```qiime2.plugin.Plugin```对象，并注册在QIIME 2框架中可发现的 *actions* ， *data formats* 和/或 *semantic types* 。  

**Result**  
*Artifact* 或 *visualizer* 的一般术语。  
*Result* 由 *method* ， *visualizer* 或 *pipeline* 产生。

**Visualization**  
可以作为QIIME *visualizer* 的输出生成的数据。  
*Visualizations* 通常在写入文件时具有扩展名```.qzv```。  

**Visualizer**  
将 *artifacts* 和 *parameters* 的某种组合作为输入的 *actions*，并且仅生成一个 *visualization* 作为输出。  
根据定义，输出 *visualizations* 不能用作其他QIIME2 *methods* 或 *visualizer* 的输入。  
*Visualizer* 只能在QIIME分析中生成终端输出。
