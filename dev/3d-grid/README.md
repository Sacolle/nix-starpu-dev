# Metas da aplicação 3d atual

- [x] Fazer uma aplicação que computa um filtro da média sobre um espaço 3d segmentado em blocos
- [x] Parametrizar o tamanho do kernel da aplicação (não necessariamente diagonal)
- [x] Obter paralelismo entre as iterações usando blocos separados
- [x] Parametrizar o número de blocos existentes
- [x] Manter a implementação com um único stencil

## Metas com a aplicação funcionando

Com a aplicação funcionando, deve-se então começar a realizar alguns experimentos. 
Por enquanto, o único tipo de kernel implementado é de CPU, tendo que expandir para uma implementação também em GPU.
Para isso é preciso de uma placa de vídeo para realizar os testes locais. Além disso, pode-se realizar experimentos no PCAD,
para avaliar o tempo de computação dado a variação dos parâmetros nesse modelo. 
Os parâmetros atuais são o **tamanho**, **segmentação do programa**, **tamanho do kernel** e **iterações**. Entretanto,
o StarPU oferece diferêntes politicas de escalonamento de programas. Dessa forma, pode-se fazer um projeto experimental
com os parâmetros acima e as possíveis políticas de escalonamento, e executar tudo isso no PCAD.
Eu preciso aprender como acessar e mexer no PCAD, mais vejo que isso virá depois ou no decorrer do processo de programar o kernel na GPU.


- [ ] Obter uma placa de vídeo.
- [ ] Implementar a computação do kernel em GPU (provavelmente CUDA).
- [ ] Montar o projeto experimenetal em R.
- [ ] Rodar os experimentos no PCAD.
