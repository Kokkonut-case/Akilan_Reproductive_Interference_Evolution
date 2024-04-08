function [data,warning]=phenologysimu_modified(popsize,d,v,C,L,mut,mutsize,initmean,shared,generations)
% simulation approach
% remember to put in the paper: address info, also the wallpaper paper

v=v/(popsize/2); % to make the situation comparable with the assumptions in the numerical stuff
warning=0; % so far no warnings found that 'trackmatings' below is too long

% initialize the population
% 'sex=1' simply defines that 1st column will sex (where we will use the notation -1 = male, +1 = female)
sex=1; pop=sign(randn([1 popsize]));
% columns 2-3 are shared alleles, 4-5 male alleles, 6-7 female alleles
sa=2:3; ma=4:5; fa=6:7; % these are indices for: shared alleles, male alleles, female alleles
maledeathrate=d(1); femaledeathrate=d(2);

trackmatings=20; % number of matings we track for a female - if she's mated with >20 males and the first 20 are all sperm depleted this is a potential problem, in that case please increase this value

pop(sa(1):fa(2),:)=randn([6 popsize]); 
pop(sa(1:2),:)=pop(sa(1:2),:)+initmean(1);
pop(ma(1:2),:)=pop(ma(1:2),:)+initmean(2);
pop(fa(1:2),:)=pop(fa(1:2),:)+initmean(3);

if shared
    mel=[sa ma]; % mel stands for male-expressed loci
    fel=[sa fa]; % fel stands for female-expressed loci
else
    mel=ma;
    fel=fa;
end

% run one mating season
% figure out when each individual emerges (according to their timing phenotype)

data=nan([((generations-1)/10)+1 4]); % we plan to plot data every 10 generations (every generation clutters the timeline with too many points...)

t_for_axis=0:10:generations;
for k=1:generations

    phenotype=(pop(sex,:)<0).*mean(pop(mel,:))+(pop(sex,:)>0).*mean(pop(fel,:));
    emergencetime=emergencesample(phenotype);

    % figure out when each individual dies (exponentially distributed time after emergence)
    deathtime=emergencetime+((pop(sex,:)<0).*exprnd(1/maledeathrate,1,popsize)+(pop(sex,:)>0).*exprnd(1/femaledeathrate,1,popsize));

    % go through all females
    f=find(pop(sex,:)>0); m=find(pop(sex,:)<0);
    T_encounter=nan*ones([trackmatings popsize]); who_encounter=nan*ones([trackmatings popsize]);
    for i=1:length(f)
        % each male's potential time of finding her
        t_encounter=max(emergencetime(f(i)),emergencetime)+exprnd(1/v,1,popsize); % both have to be alive, and have found each other
        t_encounter(pop(sex,:)>0)=NaN; % only males have encounter times with the local female
        t_encounter(t_encounter>deathtime(f(i)))=NaN; % if the female died before the encounter, the encounter never happened
        t_encounter(t_encounter>deathtime & pop(sex,:)<0)=NaN; % ditto for the male
        % store the first 'trackmatings' (usually 20) encounters. We will later check that this approximation causes
        % no error, by seeing that early encounters (lower than the 20th) usually offer sperm
        [sorted,whichmale]=sort(t_encounter);
        T_encounter(1:trackmatings,f(i))=sorted(1:trackmatings); who_encounter(1:trackmatings,f(i))=whichmale(1:trackmatings);
    end
    % check if storing trackmatings (usually 20) males was enough; issue a
    % warning if not. First flag females who mated a lot
    flag=find(~isnan(T_encounter(end,:)));

    % now figure out which matings actually happen with males who have sperm
    matecount=zeros([1 popsize]); % number of times mated (will only be used for males)
    mateidentity=NaN*ones([1 popsize]);
    relativefecundity=zeros([1 popsize]);
    while any(~isnan(T_encounter(1,:))) % nan here will denote females who have already been mated, or those who were never found by any male (who still has sperm)
        [time,ind]=min(T_encounter(1,:));
        male=who_encounter(1,ind);
        % is it possible to mate this female? yes if this male's matecount is below the allowed maximum
        if matecount(male)<C
            T_encounter(:,ind)=NaN;
            matecount(male)=matecount(male)+1;
            mateidentity(ind)=male;
            relativefecundity(ind)=400*exp(time/L)/(1+exp(time/L))^2;  % the 400 here doesn't have an impact on anything but we include it to make it match the main text of the paper
        else % we have a male, but he has no sperm
            T_encounter(:,ind)=[T_encounter(2:end,ind); NaN]; % all other mates (if they exist) are shifted upwards in the table to form the next candidates that can be checked
            who_encounter(:,ind) = [who_encounter(2:end, ind); NaN];
        end

    end

    % check if flagged (lots-mated) females actually got sperm. If not,
    % issue a warning
    if length(flag)>0 && any(relativefecundity(flag)==0) warning=1; end

    % collect data, now that matings are known
    if (k-1)/10==floor((k-1)/10)
        data((k-1)/10+1,:)=[mean(mean(pop(mel,pop(sex,:)<0))) mean(mean(pop(fel,pop(sex,:)>0))) mean(relativefecundity(relativefecundity>0)) 1-sum(~isnan(mateidentity))/sum(pop(sex,:)>0)];
        if (k-1)/100==floor((k-1)/100) 
            figure(1); 
            subplot(3,1,1); plot(t_for_axis,data(:,1),'b',t_for_axis,data(:,2),'r'); ylabel('Timing'); title(['Warning status: ' int2str(warning)])
            subplot(3,1,2); plot(t_for_axis,data(:,3),'g'); ylabel('fecundity');
            subplot(3,1,3); plot(t_for_axis,data(:,4),'k'); ylabel('matelessness'); drawnow; 
        end;
    end

    % we now have all the mating data, we can proceed to producing offspring (soft selection: their total number equals popsize)
    %
    % relative fecundity of each female depends on timing of their mating
    offspring=sign(randn([1 popsize])); % 1:1 randomly chosen sex
    offspring(sa(1):fa(2),:)=NaN;
    % mated females contribute to the next generation in proportion to
    % their fecundity (with demographic stochasticity)
    
    
    mother=randsample(popsize,popsize,true,relativefecundity);
    sire=mateidentity(mother);
    whichallele_egg=unidrnd(2,3,popsize)-1; % this randomizes between the 2 options of alleles available to transmit from mother to offspring
    whichallele_sperm=unidrnd(2,3,popsize)-1; % same for sperm
    for i=1:popsize % each female produces one clutch of offspring
        % Clarity note for below: 'male' and 'female' refer to alleles that are expressed in individuals of that sex, n
        % not to whether they are maternally or paternally inherited (i.e., either can be inherited from either parent). 
        % Whether they come from mom or dad is specified using 'sperm' or 'egg' terminology instead.
        offspring(sa(1),i)=pop(sa(1)+whichallele_egg(1,i),mother(i))';
        offspring(sa(2),i)=pop(sa(1)+whichallele_sperm(1,i),sire(i))';
        offspring(ma(1),i)=pop(ma(1)+whichallele_egg(2,i),mother(i))';
        offspring(ma(2),i)=pop(ma(1)+whichallele_sperm(2,i),sire(i))';
        offspring(fa(1),i)=pop(fa(1)+whichallele_egg(3,i),mother(i))';
        offspring(fa(2),i)=pop(fa(1)+whichallele_sperm(3,i),sire(i))';
    end
    pop(sa(1):fa(2),:)=offspring(sa(1):fa(2),:)+(rand([6 popsize])<mut).*mutsize.*randn([1 popsize]); % add mutations, and the new generation is ready
end